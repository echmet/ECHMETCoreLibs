#include "inputreader.h"

InputReader::InputReader()
{
}

InputReader::~InputReader()
{
	for (TrackedDataMap::iterator it = m_trackedData.begin(); it != m_trackedData.end(); it++) {
		freeData(&it->second);
	}
}

void InputReader::freeData(const constituent_array_t *array)
{
	ldr_destroy_constituent_array(array);
}

constituent_array_t InputReader::read(const std::string &filepath)
{
	enum LoaderErrorCode errorCode = JLDR_OK;
	const constituent_array_t ctuents = ldr_load_constituents(filepath.c_str(), &errorCode);

	if (ctuents.constituents == NULL) {
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
	m_trackedData[filepath] = ctuents;

	return ctuents;
}

void InputReader::release(const std::string &filepath)
{
	TrackedDataMap::iterator it = m_trackedData.find(filepath);

	if (it != m_trackedData.end()) {
		freeData(&it->second);
		m_trackedData.erase(it);
	}
}
