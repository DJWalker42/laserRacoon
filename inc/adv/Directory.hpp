#ifndef DIR_INCLUDED
#define DIR_INCLUDED

#include <list>
#include <string>

namespace dir{

	std::string current_dir();
		// returns the current working directory as a string

	bool change_dir(const std::string &new_dir);
		// returns true if successfully changed to new_dir

	void list_dir(std::list<std::string> &files, const std::string &filespec = "*.*");
		// fills files with listing of current directory

	bool make_dir(const std::string &new_dir);
		// returns true if successfully creates new_dir

	bool del_dir(const std::string &old_dir);
		// returns true if successfully deletes old_dir

}

#endif
