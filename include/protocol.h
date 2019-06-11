#pragma once
#include <sstream>
#include <string>
#include <vector>

namespace protocol{

pendulum_on_cart deserialize(const std::string &str){
	std::vector<double> state;
	std::vector<std::string> strings;
	{
		std::stringstream ss(str);
		std::string buffer;

		while(std::getline(ss, buffer, ',')){
			strings.emplace_back(buffer);
		}
	}

	id = std::stoi(strings.at(0));
	w = std::stod(strings.at(1));
	l = std::stod(strings.at(2));
	state.reserve(strings.size()-3);
	for(int i = 3; i<strings.size(); i++){
		state.emplace_back(std::stod(strings.at(i)));
	}
	return pendulum_on_cart(state);
}
};
