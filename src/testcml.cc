
#include <iostream>
#include "xml.h++"

void printXMLList(std::list<XML> &);
std::list<XML> parseXMLfile(char * filename);

int main(int argc,char *argv[])
{
  std::list<XML> lxml;
  lxml = parseXMLfile(argv[1]);
  printXMLList(lxml);

  return 0;
}
