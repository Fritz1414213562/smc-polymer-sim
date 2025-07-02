#include "Utility.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>


namespace Utility
{

std::string get_file_suffix(const std::string& filename)
{
	const std::size_t file_suffix_from = filename.rfind(".");
	if (file_suffix_from == std::string::npos)
		throw std::runtime_error(
			"[error] There is no file extension in " + filename + "."
			" The file type can not be specified");
	const std::size_t file_suffix_len = filename.length() - file_suffix_from;
	return filename.substr(file_suffix_from, file_suffix_len);
}

std::string erase_space(std::string&& str)
{
	const auto new_end = std::remove_if(str.begin(), str.end(),
							[](const char x){ return std::isspace(x); });
	str.erase(new_end, str.end());
	return str;
}

void clear_file(const std::string& filename)
{
	std::filesystem::path fpath(filename);
	fpath.remove_filename();
	if (fpath.empty())
	{
		fpath = "./";
	}

	if (! std::filesystem::exists(fpath))
		std::filesystem::create_directories(fpath);
	if (! std::filesystem::exists(fpath) || ! std::filesystem::is_directory(fpath))
		throw std::runtime_error(
			"failed to make output directory: " + fpath.string());
	std::ofstream ofs(filename);
	if (not ofs.good())
		throw std::runtime_error("file open error : " + filename);
	ofs.close();
	return ;

}

const toml::value& find_either(const toml::value& v,
        const std::string& key1, const std::string& key2)
{
    // A functor to find a value that corresponds to either of the key.
    // If both key exists, throw an error.
    if(v.contains(key1) && v.contains(key2) != 0)
    {
        std::cerr << toml::format_error("[error] key duplicates.", v.at(key1), "here", v.at(key2),
                                        "this conflicts with the above value definition")
                 << std::endl;
    }
    if(v.contains(key1)) { return v.at(key1); }
    if(v.contains(key2)) { return v.at(key2); }

    throw std::runtime_error(toml::format_error("both keys, \"" + key1 + "\" and \"" + key2 +
                             "\", are not found.", v, "in this table"));
}

void merge_toml_tables(toml::value& table, const toml::value& other)
{
    assert(table.is_table());
    assert(other.is_table());

    for(const auto& kv : other.as_table())
    {
        if(table.contains(kv.first))
        {
            if(table.at(kv.first).is_table())
            {
                merge_toml_tables(table.at(kv.first), kv.second);
            }
            else if(table.at(kv.first).is_array() && kv.second.is_array())
            {
                toml::array        arr       = table.at(kv.first).as_array();
                const toml::array& other_arr = kv.second.as_array();
                arr.insert(arr.end(), other_arr.begin(), other_arr.end());
                table.at(kv.first) = arr;
            }
            else
            {
                throw std::runtime_error(toml::format_error("value \"" +
                    kv.first + "\" duplicates.", table, "first defined here",
                    kv.second, "and also defined here"));
            }
        }
        else
        {
            table.as_table().emplace(kv.first, kv.second);
        }
    }
}

void expand_include(toml::value& v)
{
    if(!v.is_table() && !v.is_array()) {return;}

    if(v.is_table())
    {
        for(auto& kv : v.as_table())
        {
            expand_include(kv.second);
        }

        // expand include in this table
        if(v.contains("include"))
        {
            if(v.at("include").is_array())
            {
                for(auto fname : toml::find<std::vector<std::string>>(v, "include"))
                {
                    std::cerr << "    expanding file " << fname << std::endl;
                    merge_toml_tables(v, toml::parse(fname));
                }
            }
            else
            {
                const auto& fname = toml::find<std::string>(v, "include");
                std::cerr << "    expanding file " << fname << std::endl;
                merge_toml_tables(v, toml::parse(fname));
            }
            v.as_table().erase("include");
        }
    }
    else if(v.is_array()) // handle an array of tables
    {
        for(auto& elem : v.as_array())
        {
            expand_include(elem);
        }
    }
    return;
}

void add_offset(std::size_t& index, const toml::value& offset)
{
    if(offset.is_array())
    {
        throw std::runtime_error("[error] add_offset: "
            "offset array is invalid here. expected type is scalar.");
    }
    else
    {
        index += offset.as_integer();
    }
}

template <typename TupleLike>
void add_offset_impl(TupleLike& indices, const toml::value& offset)
{
    if(offset.is_array())
    {
        if(offset.size() != std::tuple_size_v<TupleLike>)
        {
            throw std::runtime_error("[error] add_offset: "
                "invalid size of offset array " + std::to_string(offset.size()) +
                 " here. the size of offset array must much to the size of "
                 "indices array. expected size is " +
                 std::to_string(std::tuple_size_v<TupleLike>) + ".");
        }
        Utility::for_each(indices, [&offset](const std::size_t i, auto& elem) {
                elem += offset.at(i).as_integer();
            });
    }
    else
    {
        Utility::for_each(indices, [&offset](const std::size_t i, auto& elem) {
                elem += offset.as_integer();
            });
    }
}

void add_offset(std::pair<std::size_t, std::size_t>& indices, const toml::value& offset)
{
    add_offset_impl(indices, offset);
}
void add_offset(std::array<std::size_t, 3>& indices, const toml::value& offset)
{
    add_offset_impl(indices, offset);
}
void add_offset(std::array<std::size_t, 4>& indices, const toml::value& offset)
{
    add_offset_impl(indices, offset);
}

bool check_keys_available(const toml::value& table,
                                 std::initializer_list<std::string> list)
{
    bool all_available = true;
    for(const auto& kv : table.as_table())
    {
        if(list.end() == std::find(list.begin(), list.end(), kv.first))
        {
            std::cerr << "\033[33m[warning]\033[m"
                      << " unknown value \"" << kv.first << "\" found. this "
                      << kv.second.type() << " will never be used."
                      << std::endl;
            all_available = false;
        }
    }

    return all_available;
}

}
