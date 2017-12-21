#ifndef ADDITIONAL_TOOLS_H
#define ADDITIONAL_TOOLS_H
#include <vector>
#include <iostream>
#include <map>

class SymbolChange;

template<class key,class value>
std::ostream& operator <<(std::ostream &out, const std::map<key,value> &map );

std::ostream& operator << (std::ostream &out, const SymbolChange &symb );

#include "additional_tools.tpp"
#endif