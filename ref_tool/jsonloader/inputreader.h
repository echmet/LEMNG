#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <map>
#include <string>
#include "constituents_json_ldr.h"

//typedef struct constituent_array constituent_array_t;

class InputReader {
public:
	class InputReaderException : public std::exception {
	public:
		const char * what() const throw()
		{
			return "Unspecified input reader exception";
		}
	};

	class FileErrorException : public InputReaderException {
	public:
		const char * what() const throw()
		{
			return "Unable to open the file for reading";
		}
	};

	class InvalidInputException : public InputReaderException {
	public:
		const char * what() const throw()
		{
			return "System composition is invalid";
		}
	};

	class MalformedInputException : public InputReaderException {
	public:
		const char * what() const throw()
		{
			return "Input file is malformed";
		}
	};

	class NoMemoryException : public InputReaderException {
	public:
		const char * what() const throw()
		{
			return "Insufficient memory";
		}
	};

	InputReader();
	~InputReader();
	const constituent_array_t * read(const std::string &filepath);
	void release(const std::string &filename);

private:
	typedef std::map<std::string, const constituent_array_t *> TrackedDataMap;

	void freeData(const constituent_array_t *array);

	TrackedDataMap m_trackedData;
};

#endif // INPUTREADER_H
