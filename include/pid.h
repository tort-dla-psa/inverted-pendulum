#pragma once
#include <cmath>

class pid{
	double kp, ki, kd;
	double pgain, igain, dgain;
	double prev_actual, prev_output;
	double error_sum;
	double min_output, max_output;

	inline void norm(double &var, double min, double max){
		if(var < min){
			var = min;
		}else if(var > max){
			var = max;
		}
	}
public:
	pid(double kp, double ki, double kd, double min_output, double max_output) {
		this->min_output = min_output;
		this->max_output = max_output;
		this->kp = kp;
		this->ki = ki;
		this->kd = kd;
		error_sum = 0;
	}

	void set_kp(double kp){ this->kp = kp; }
	void set_ki(double ki){ this->ki = ki; }
	void set_kd(double kd){ this->kd = kd; }
	double get_kp()const{ return kp; }
	double get_ki()const{ return ki; }
	double get_kd()const{ return kd; }
	double get_summ()const{ return std::abs(error_sum); }
	double get_min_force()const{ return min_output; }
	double get_max_force()const{ return max_output; }

	double update(double target, double current){
		const double error = target - current;

		pgain = kp * error;
		if(error_sum == 0){
			prev_actual = current;
			prev_output = pgain;
		}

		dgain = -kd*(current - prev_actual);
		prev_actual = current;

		igain = ki * error_sum;
		error_sum += error;

		double output = pgain + igain + dgain;
		if(output > max_output || output < min_output){
			norm(output, min_output, max_output);
		}
		prev_output = output;
		return output;
	}
};
