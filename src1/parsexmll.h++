/*! \file parsexmll.h++
  holds the functions for parsing a xml files
*/


#ifndef _PARSEXMLL_H
#define _PARSEXMLL_H

std::list<XML> parseXMLfile(const std::string &filename);
// std::list<XML> clist(std::list<XML>::iterator p);
void printXMLList(std::list<XML>::iterator it,std::list<XML>::iterator end);
void printXMLListl(std::list<XML>::iterator it,std::list<XML>::iterator end);
void printXMLList(std::list<XML> &);
#endif
