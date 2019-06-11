#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <random>
#include <vector>
#include <memory>
#include <atomic>
#include <sys/socket.h>
#include <boost/numeric/odeint.hpp>
#include "socket_op.h"
#include "file.h"
#include "pid.h"
#include <pugixml.hpp>

#define pi 3.1415926
#define g 9.81

typedef std::vector<double> state_vec;

//x, x_speed, rad, rad_speed
state_vec init = {0., 0., 1., 0.};

int x_place = 0,
	x_place_speed = 1,
	angle_place = 2,
	angle_speed_place = 3;

class pendulum_on_cart{
protected:
	double l, m_pend, m_cart;
	double force;
	state_vec state;
	double x_speed_decr, rad_speed_decr, force_decr;
public:
	pendulum_on_cart(double l, double m_pend, double m_cart, const state_vec &state,
		double x_speed_decr, double rad_speed_decr, double force_decr)
	{
		force = 0.0;
		this->l = l;
		this->m_pend = m_pend;
		this->m_cart = m_cart;
		this->state = state;
		this->x_speed_decr = x_speed_decr;
		this->force_decr = force_decr;
		this->rad_speed_decr = rad_speed_decr;
	}

	pendulum_on_cart()
		:pendulum_on_cart(1, 1, 1, {0., 0., 0., 0.}, 0, 0, 0)
	{}

	~pendulum_on_cart(){}

	void operator() (const state_vec x, state_vec &dxdt, const double){
		double _sin = std::sin(x[2]), _cos = std::cos(x[2]);
		dxdt[0] = x[1];
		dxdt[1] = force/m_cart - l*x[3]*x[3]*_cos + g*_cos*_sin;
		dxdt[1] /= (-_sin*_sin + m_cart/m_pend + 1.);
		dxdt[2] = x[3];
		dxdt[3] = (-g - 1.*_sin*_cos*dxdt[1])/l*_cos;

		dxdt[1] -= x_speed_decr*x[1];
		dxdt[3] -= rad_speed_decr*x[3];
	}

	double get_l()const{ return l; }
	double get_rad()const{ return state[angle_place]; }
	double get_rad_speed()const{ return state[angle_speed_place]; }
	double get_x()const{ return state[x_place]; }
	double get_x_speed()const{ return state[x_place_speed]; }
	double get_m_cart()const{ return m_cart; }
	double get_m_pend()const{ return m_pend; }
	double get_force()const{ return force; }
	state_vec get_state()const{ return state; }

	void set_force(double force){ this->force = force; }
	void set_state(const state_vec &state){ this->state = state; }
};

template <typename T>
class simulation{
protected:
	T obj;
	pid p;
	boost::numeric::odeint::runge_kutta_dopri5<state_vec> stepper;

	std::thread* sim_thr = nullptr;
	std::chrono::milliseconds sim_sleep;
	double step;
	std::atomic_bool end_requested;

	std::thread* pid_thr = nullptr;
	std::mutex pid_mt;
	int pid_hz;

	void sim_func(){
		double i=0;
		while(!end_requested){
			pid_mt.lock();
			auto state = obj.get_state();
			stepper.do_step(obj, state, i, step);
			auto& ref = state[angle_place];
			while(ref > pi){
				ref -= pi*2.;
			}
			while(ref < -pi){
				ref += pi*2.;
			}
			obj.set_state(state);
			pid_mt.unlock();
			i += step;
			std::this_thread::sleep_for(sim_sleep);
		}
	}

	void pid_func(){
		while(!end_requested){
			pid_mt.lock();
			double out = p.update(pi/2., obj.get_rad());
			obj.set_force(out);
			pid_mt.unlock();
			auto now = std::chrono::steady_clock::now();
			auto next = now + std::chrono::milliseconds((int)(1000./pid_hz));
			std::this_thread::sleep_until(next);
		}
	}

public:
	simulation(T obj, pid p, std::chrono::milliseconds sim_sleep, double step, int pid_hz)
		:p(p),
		obj(obj),
		sim_sleep(sim_sleep),
		step(step),
		pid_hz(pid_hz)
	{}

	~simulation(){
		stop();
		delete pid_thr;
		delete sim_thr;
	}

	pid get_pid()const{
		return p;
	}

	void start(){
		end_requested = false;
		delete pid_thr;
		delete sim_thr;
		pid_thr = new std::thread(&simulation::pid_func, this);
		sim_thr = new std::thread(&simulation::sim_func, this);
	}

	void stop(){
		end_requested = true;
		if(sim_thr && sim_thr->joinable()){
			sim_thr->join();
		}
		if(pid_thr && pid_thr->joinable()){
			pid_thr->join();
		}
	}

	state_vec get_state(){
		auto bak = obj.get_state();
		return bak;
	}
};

template<typename T>
class multiple_sim{
	std::vector<std::shared_ptr<simulation<T>>> sims;
	typedef std::vector<std::tuple<double, double, double>> coefs_vec;

	static coefs_vec make_coefs(size_t count, double kp, double ki, double kd,
		double percent)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(-percent, percent);

		coefs_vec vec;
		vec.reserve(count);
		double newkp = kp;
		double newki = ki;
		double newkd = kd;
		vec.emplace_back(std::make_tuple(newkp, newki, newkd));
		for(size_t i=1; i<count; i++){
			double newkp = kp + kp*dis(gen);
			double newki = ki + ki*dis(gen);
			double newkd = kd + kd*dis(gen);
			vec.emplace_back(std::make_tuple(newkp, newki, newkd));
		}
		return vec;
	}
public:
	multiple_sim(T prototype, pid p, size_t count,
		std::chrono::milliseconds sim_sleep, double sim_step,
		size_t pid_hz, double percent)
		:multiple_sim(prototype, sim_sleep, sim_step, pid_hz,
			make_coefs(count, p.get_kp(), p.get_ki(), p.get_kd(), percent),
			p.get_min_force(), p.get_max_force())
	{}

	multiple_sim(T prototype, std::chrono::milliseconds sim_sleep, double sim_step,
		size_t pid_hz, const coefs_vec &coef, double min_force, double max_force)
	{
		const size_t size = coef.size();
		sims.reserve(size);

		for(auto &c:coef){
			pid p(std::get<0>(c), std::get<1>(c), std::get<2>(c),
				min_force, max_force);
			sims.emplace_back(new simulation<T>(prototype, p,
				sim_sleep, sim_step, pid_hz));
		}
	}

	std::vector<std::shared_ptr<simulation<T>>> get_sims()const{
		return sims;
	}
	
	void start(){
		for(auto &sim:sims){
			sim->start();
		}
	}

	void stop(){
		for(auto &sim:sims){
			sim->stop();
		}
	}
};


void print(const std::unique_ptr<class pid> &pid){
	std::cout<<"p"<<pid->get_kp()<<"|"
		<<"i"<<pid->get_ki()<<"|"
		<<"d"<<pid->get_kd()<<"|"
		<<"summ"<<pid->get_summ()<<"\n";
}

std::shared_ptr<multiple_sim<pendulum_on_cart>> sim;
std::atomic_bool sim_run;

template<typename T>
void sim_func(T time){
	sim->start();
	std::this_thread::sleep_for(time);
	sim->stop();
	sim_run = false;
}

void check(pugi::xml_parse_result result, const std::string &source){
	if(!result){
		std::cout << "Error description: " << result.description() << "\n";
		std::cout << "Error offset: " << result.offset
			<< " (error at [..." << (source + std::to_string(result.offset))
			<< "]\n\n";
	}
}

void print(pugi::xml_node node){
	std::cout<<"node name: "<<node.name()<<"\n";
	std::cout<<"node children: ";
	auto ch = node.children();
	for(auto c:ch){
		std::cout<<c.name()<<" ";
	}
	std::cout<<"\n";
}

template<typename T>
void fill(pugi::xml_node node, T &val){
	if(std::is_same<T, double>())
		val = std::stod(node.child_value());
	else
		val = std::stoi(node.child_value());
	if(node){
		std::cout<<"filling from node: "<<node.name()
			<<" val: "<<node.child_value()<<"\n";
	}else{
		std::cerr<<"filling from null node!"<<"\n";
	}
}

int main(){
	const std::string path = "/tmp/sim";
	const auto conf_path = "conf.xml";
	double cart_w, pole_l, m_cart, m_pole, step=0.1;
	double rad_speed_decr, x_speed_decr, force_decr;
	std::chrono::milliseconds sim_sleep, sleep_time, sample_time;
	double kp, ki, kd, force_max, force_min;
	double init_perc, perc;
	size_t size;
	size_t pid_hz;

	{
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(conf_path);
		check(result, conf_path);
		auto root = doc.document_element();
		print(root);
		fill(root.child("count"), size);
		fill(root.child("pid_hz"),pid_hz);
		fill(root.child("cart").child("width"),cart_w);
		fill(root.child("cart").child("mass"),m_cart);
		fill(root.child("pole").child("length"),pole_l);
		fill(root.child("pole").child("mass"),m_pole);
		{
			auto ode = root.child("ode");
			print(ode);
			fill(ode.child("rad_speed_decr"), rad_speed_decr);
			fill(ode.child("x_speed_decr"), x_speed_decr);
			fill(ode.child("force_decr"), force_decr);
		}
		{
			auto sleeps = root.child("sleeps");
			print(sleeps);
			int temp;
			fill(sleeps.child("phys_simulator"), temp);
			sim_sleep = std::chrono::milliseconds(temp);
			fill(sleeps.child("gen_algorithm"), temp);
			sleep_time = std::chrono::milliseconds(temp);
			fill(sleeps.child("sample"), temp);
			sample_time = std::chrono::milliseconds(temp);
		}
		{
			auto pid = root.child("pid");
			print(pid);
			fill(pid.child("kp"), kp);
			fill(pid.child("ki"), ki);
			fill(pid.child("kd"), kd);
			fill(pid.child("force_min"), force_min);
			fill(pid.child("force_max"), force_max);
		}
		{
			auto percents = root.child("percents");
			print(percents);
			fill(percents.child("init_pid_perc"), init_perc);
			fill(percents.child("pid_perc"), perc);
		}
	}
	

	IO::socket_op s_op;
	auto sock = s_op.create(path, AF_UNIX, SOCK_STREAM, 0);
	s_op.bind(sock);
	s_op.listen(sock, 1);
	auto cli = s_op.accept(sock);

	pendulum_on_cart pend(pole_l, m_cart, m_pole, init,
		x_speed_decr, rad_speed_decr, force_decr);
	std::shared_ptr<pid> best = nullptr;

	while(true){
		if(!best){
			pid temp(kp, ki, kd, force_min, force_max);
			sim = std::make_shared<multiple_sim<pendulum_on_cart>>
				(pend, temp, size, sim_sleep, step, pid_hz, init_perc);
		}else{
			sim = std::make_shared<multiple_sim<pendulum_on_cart>>
				(pend, *best, size, sim_sleep, step, pid_hz, perc);
		}
		sim_run = true;
		std::thread sim_thr(&sim_func<std::chrono::milliseconds>, sleep_time);
		
		while(sim_run){
			auto sims = sim->get_sims();
			const size_t size = sims.size();

			s_op.send(cli, (int)1);
			s_op.send(cli, "c");
			for(size_t i=0; i<size; i++){
				std::string str = "d" + std::to_string(i) + ",";
				auto state = sims[i]->get_state();
				str += std::to_string(cart_w) + ",";
				str += std::to_string(pole_l) + ",";
				state[angle_place] *= -1;
				for(const auto &var:state){
					str += std::to_string(var) + ",";
				}
				//std::cout<<str<<'\n';
				s_op.send(cli, (int)str.size());
				s_op.send(cli, str);
			}
			std::this_thread::sleep_for(sample_time);
		}
		sim_thr.join();
		auto sims = sim->get_sims();

		const size_t size = sims.size();

		if(!best){
			best = std::make_shared<pid>(sims[0]->get_pid());
		}
		bool best_set = sims[0]->get_state()[angle_place] > 0;
		for(const auto &sim:sims){
			auto newp = sim->get_pid();
			auto state = sim->get_state();
			if(newp.get_summ() < best->get_summ()
				&& state[angle_place] > 0)
			{
				best = std::make_shared<pid>(newp);
				best_set = true;
			}
		}
		if(!best_set){
			best = nullptr;
			std::cout <<"best selection failed\n";
		}else{
			std::cout <<"best pid kp:"<<best->get_kp()
				<<" ki:"<<best->get_ki()
				<<" kd:"<<best->get_kd()
				<<" error:"<<best->get_summ()<<"\n";
		}
	}
	try{
		s_op.close(cli);
	}catch(...){};
	s_op.close(sock);
}
