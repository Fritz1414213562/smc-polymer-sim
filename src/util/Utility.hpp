#ifndef POLYMER_MDMC_UTILITY_HPP
#define POLYMER_MDMC_UTILITY_HPP

#include <toml.hpp>
#include <iostream>
#include <string>
#include <toml_fwd.hpp>

namespace Utility
{

// file output
template<typename T>
void write_as_bytes(std::ofstream& ofs, const T& v) noexcept
{
	using Type = typename std::remove_reference<T>::type;
	ofs.write(reinterpret_cast<const char*>(std::addressof(v)), sizeof(Type));
	return;
}

std::string get_file_suffix(const std::string& filename);
std::string erase_space(std::string&& str);
void clear_file(const std::string& filename);

// file input
void merge_toml_tables(toml::value& table, const toml::value& other);
void expand_include(toml::value& v);
void add_offset(std::size_t& index, const toml::value& offset);
void add_offset(std::pair<std::size_t, std::size_t>& indices, const toml::value& offset);

void add_offset(std::array<std::size_t, 3>& indices, const toml::value& offset);
void add_offset(std::array<std::size_t, 4>& indices, const toml::value& offset);

// This check all the keys in a table are found in a list.
//     If there is a key that is not found in the range, it warns about the
// corresponding value will be ignored.
//
// Use it as the following.
// ```cpp
// check_keys_available(table, {"foo", "bar", "baz"});
// ```
bool check_keys_available(const toml::value& table,
                          std::initializer_list<std::string> list);

}


#endif // POLYMER_MDMC_UTILITY_HPP
