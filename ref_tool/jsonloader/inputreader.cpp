#include "inputreader.h"

InputReader::InputReader()
{
}

InputReader::~InputReader()
{
	for (TrackedDataMap::iterator it = m_trackedData.begin(); it != m_trackedData.end(); it++) {
		freeData(it->second);
	}
}

void InputReader::freeData(const constituent_array_t *array)
{
	for (size_t idx = 0; idx < array->count; idx++)
		ldr_destroy_constituent(&array->constituents[idx]);

	ldr_destroy_array(array);
}

const constituent_array_t * InputReader::read(const std::string &filepath)
{
	enum LoaderErrorCode errorCode = JLDR_OK;
	const constituent_array_t *array = ldr_loadFromFile(filepath.c_str(), &errorCode);

	if (array == NULL) {
		switch (errorCode) {
		case JLDR_E_BAD_INPUT:
			throw InvalidInputException();
			break;
		case JLDR_E_CANT_READ:
			throw FileErrorException();
			break;
		case JLDR_E_MALFORMED:
			throw MalformedInputException();
			break;
		case JLDR_E_NO_MEM:
			throw NoMemoryException();
			break;
		default:
			throw InputReaderException();
			break;
		}
	}

	release(filepath);
	m_trackedData[filepath] = array;

	return array;
}

void InputReader::release(const std::string &filepath)
{
	TrackedDataMap::iterator it = m_trackedData.find(filepath);

	if (it != m_trackedData.end()) {
		freeData(it->second);
		m_trackedData.erase(it);
	}
}
