//
// Created by gharg on 25.09.17.
//

#ifndef PROJECT_TESTINSTANCECONTAINER_H
#define PROJECT_TESTINSTANCECONTAINER_H

#include <iostream>
#include <fstream>
#include "string"
#include "boost/filesystem.hpp"
#include <boost/filesystem/fstream.hpp>
#include "TestInstance.hpp"

namespace MacroTesting {
	class TestInstanceContainer {
		std::vector<TestInstance<NT>> instances;
	public:
		void importInstances(std::string path) {
			std::cout << "Path: " << path << std::endl;
			boost::filesystem::path p(path);
			if (boost::filesystem::exists(p)) {
				std::cout << "file/directory exists" << std::endl;
				if (boost::filesystem::is_directory(p)) {
					for (boost::filesystem::directory_entry &it : boost::filesystem::directory_iterator(p)) {
						instances.push_back(readInstance(it.path()));
					}

				} else if (boost::filesystem::is_regular_file(p)) {

					instances.push_back(readInstance(p));
				}
			}
		}

		std::vector<TestInstance<NT> *> getTestInstances(std::string namePattern) {
			std::vector<TestInstance<NT> *> result;
			std::string name = namePattern;
			std::transform(name.begin(), name.end(), name.begin(), ::tolower);
			//TODO: Real name pattern matching!
			std::cout << "overall instances: " << instances.size() << std::endl;
			for (uint i = 0; i < instances.size(); ++i) {
				std::string instanceName = instances[i].name;
				std::transform(instanceName.begin(), instanceName.end(), instanceName.begin(), ::tolower);
				if (instanceName.find(name) != std::string::npos) {
					result.push_back(&instances[i]);
				}
			}
			return result;
		}

		std::vector<TestInstance<NT> *> getTestInstances(ulong rangeBegin, ulong rangeEnd, ulong maxAmount = 0) {
			std::vector<TestInstance<NT> *> result;
			ulong amountCounter = ulong() - 1;
			if(maxAmount != 0) {
				amountCounter = maxAmount;
			}
			for (uint i = 0; (i < instances.size() && amountCounter > 0); ++i) {
				if (instances[i].get_point_amount() >= rangeBegin && instances[i].get_point_amount() <= rangeEnd) {
					result.push_back(&instances[i]);
					--amountCounter;
				}
			}
			return result;
		};
	private:
		TestInstance<NT> readInstance(boost::filesystem::path filename) {
			const char *lineseparator = "\n";
			const char *columnseparator = "\t";
			std::vector<NT> points;
			std::ifstream ifstream(filename.c_str());
			if (ifstream) {
				// get length of file:
				ifstream.seekg(0, ifstream.end);
				unsigned long length;
				length = (unsigned long) ifstream.tellg();
				ifstream.seekg(0, ifstream.beg);

				char *buffer = new char[length];

				//std::cout << "Reading " << length << " characters... " << std::endl;
				// read data as a block:
				ifstream.read(buffer, length);

				//if (ifstream)
				//	std::cout << "all characters read successfully." << std::endl;
				//else
				//	std::cout << "error: only " << ifstream.gcount() << " could be read" << std::endl;
				ifstream.close();

				// ...buffer contains the entire file...

				std::string csvString(buffer, length);
				// split string into substrings of lines
				unsigned long pos = 0;
				unsigned long endPos = csvString.find(lineseparator, pos);
				if (endPos == std::string::npos) {
					endPos = csvString.size() - 1;
				}
				while (pos < csvString.size()) {
					// ignore comment lines
					if (csvString[pos] != '#') {
						unsigned long innerPosMid = csvString.find(columnseparator, pos);
						std::string subStrX = csvString.substr(pos, innerPosMid);
						std::string subStrY = csvString.substr(innerPosMid + 1, endPos);
						double x = std::stod(subStrX);
						double y = std::stod(subStrY);
						points.push_back(x);
						points.push_back(y);
					}


					pos = endPos + 1;
					endPos = csvString.find(lineseparator, pos);
					if (endPos == std::string::npos) {
						endPos = csvString.size() - 1;
					}
				}
				delete[] buffer;
			} else {
				std::cout << "file could not be read" << std::endl;
			}
			std::string name = filename.filename().stem().string();
			return TestInstance<NT>(name, &points[0], &points[points.size()], 2);
		}
	};
}
#endif //PROJECT_TESTINSTANCECONTAINER_H
