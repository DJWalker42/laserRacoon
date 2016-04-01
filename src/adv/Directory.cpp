#include <direct.h>
#include <io.h>
#include <cstdio>
#include <adv/Directory.hpp>

namespace dir {
	std::string current_dir()
	{
		char buffer[FILENAME_MAX] = "";
		_getcwd(buffer, sizeof buffer);
		return buffer;
	}

	bool change_dir(const std::string &new_dir)
	{
		return _chdir(new_dir.c_str()) == 0;
	}

	void list_dir(std::list<std::string> &files, const std::string &filespec)
	{
		files.clear();

		_finddata_t info;
		long cursor = _findfirst(filespec.c_str(), &info);

		if (cursor != -1)
		{
			do
			{
				files.push_back(info.name);
			} while (_findnext(cursor, &info) == 0);

			_findclose(cursor);
		}
	}

	bool make_dir(const std::string &new_dir)
	{
		return _mkdir(new_dir.c_str()) == 0;
	}

	bool del_dir(const std::string &old_dir)
	{
		return _rmdir(old_dir.c_str()) == 0;
	}

}