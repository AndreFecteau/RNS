#ifndef SERIALIZE_H
#define SERIALIZE_H

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
#include <fstream>

template<typename T>
void serialize_to_file(T& solver, std::string filename)
{
  std::ofstream os(filename + ".cereal", std::ios::binary);
  cereal::BinaryOutputArchive archive_o( os );
  archive_o(solver);
}

template<typename T>
void unserialize_to_file(T& solver, std::string filename)
{
  std::ifstream is(filename + ".cereal", std::ios::binary);
  cereal::BinaryInputArchive archive_i( is );
  archive_i(solver);
}

#endif //#ifndef SERIALIZE_H
