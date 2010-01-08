/*! 
  \file parsexmll.c++ 
  subroutine to parse a file an xml file and 
  return a tree with the parsed data
*/

#include "xml.h++"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stack>
#include <list>

const char * HEAD = "<?xml version=\"1.0\"?>";
// using namespace std;

// return the childs lists
// std::list<XML> clist(std::list<XML>::iterator p)
// {
//   if((*p).xmllist().begin() == (*p).xmllist().end()){
//     std::cerr << "NO list!\n";
//     exit(1);
//   }
//   return (*p).xmllist();
// }

// print out the list (not the sublists)
void printXMLListl(std::list<XML>::iterator it,std::list<XML>::iterator end)
{
  while(it!=end) {
    std::cout<<"<" <<(*it).tag()<<" att=\""<<(*it).att()<<"\">\n";
    {
      int j=-1,i=0;
      std::list<std::string>::iterator lnow,lend;
      lnow = (*it).databegin(); lend = (*it).dataend(); 
      std::list<XML>::iterator lxml,elxml;
      lxml = (*it).xmllist().begin();
      elxml = (*it).xmllist().end();
      for(;lnow != lend; ++lnow,++i){
	std::cout <<(*lnow)<< std::endl;
	if(lnow == (*it).pos()) j=i;
	// if(lxml!= elxml) printXMLList(lxml,++lxml);
      }
      // std::cout <<"\tpos at = "<< j << std::endl;
    }
    std::cout<<"</"<<(*it).tag()<<">\n";
    ++it;
  }
}

// print out the list (and sublists)
void printXMLList(std::list<XML>::iterator it,std::list<XML>::iterator end)
{
  while(it!=end) {
    std::cout<<"<" <<(*it).tag()<<" att=\""<<(*it).att()<<"\">";
    {
      int j=-1,i=0;
      std::list<std::string>::iterator lnow,lend;
      lnow = (*it).databegin(); lend = (*it).dataend(); 
      std::list<XML>::iterator lxml,elxml;
      lxml = (*it).xmllist().begin();
      elxml = (*it).xmllist().end();
      for(;lnow != lend; ++lnow,++i){
	std::cout <<(*lnow); // << std::endl;
	if(lnow == (*it).pos()) j=i;
	if(lxml!= elxml) {std::cout << std::endl;printXMLList(lxml,++lxml);}
      }
      // std::cout <<"\tpos at = "<< j << std::endl;
    }
    std::cout<<"</"<<(*it).tag()<<">\n";
    ++it;
  }
}

// print out the list (and sublists)
void printXMLList(std::list<XML> &src)
{
  std::list<XML>::iterator it = src.begin();
  std::cout<<"head = " << (*it).tag() << std::endl;
  std::list<std::string>::iterator lnow,lend;
  lnow = (*it).databegin(); lend = (*it).dataend(); 
  for(;lnow != lend; ++lnow) std::cout <<(*lnow)<< std::endl;
  std::list<XML>::iterator lxml,elxml;
  lxml = (*it).xmllist().begin();
  elxml = (*it).xmllist().end();
  printXMLList((*it).xmllist().begin(),(*it).xmllist().end());
  std::cout<<"\nEnd XML\n\n";
}

// skip all XML comments <!-- .... -->
static std::string skipXMLcomment(std::ifstream &file, std::string &input)
{
  std::string::size_type i1 = input.find("<");
  if(i1 == std::string::npos) return input;

  std::string sstart;
  if(i1!=0) sstart = input.substr(0,i1-1);
  std::string::size_type i2;
  
  i1 = input.find("<!--"); // is it a comment ?
  while(i1 != std::string::npos){ // we have a comment
    i2 = input.find("-->"); // end comment ?
    while(i2 == std::string::npos){
      if(getline(file,input)){
	i2 = input.find("-->"); // end comment ?	    
      } else {
	std::cerr << "Comment never finished!\n";
	exit(1);
      }
    }
    // std::cout << "found comment :" << input << std::endl;
    input.erase(0,i2+3); // continue after comment
    i1 = input.find("<!--"); // more comments ?
  } // end comment
  return sstart + input;
}

// read in an XML start tag
static int startXMLtag(std::string &sline,std::string &stag,std::string &satt)
{
  int iend=0;
  std::string::size_type i1 = sline.find("<");
  if(i1 == std::string::npos){std::cerr<<"Error! no < \n";exit(1);}
  std::string::size_type i2 = sline.find("/>");
  if(i2 != std::string::npos){ // we are at the end of the key ie < key />
    sline.erase(i2,1);
    iend=1;
  } 

  i2 = sline.find(">");
  if(i2 == std::string::npos){std::cerr<<"Error! < but no end >\n";exit(1);}
  std::string sdata = sline.substr(i1+1,i2-i1-1);
  sline.erase(0,i2-i1+1);
  std::istringstream ist(sdata);
  ist >> stag;
  // std::cout << "key = " << stag << std::endl;
  satt.erase();
  getline(ist,satt); 
  return iend;
}

// we found an error in the end tag
static void errendtag(std::string & stag,
		      std::list<XML>::iterator & trit,
		      const std::list<XML>::iterator & top)
{
  if(stag != (*trit).tag()){ // make sure we match
    std::cerr << "end tag \""<<stag<<"\" doesn't match start tag \"";
    std::cerr << (*trit).tag() << "\"\n";
    exit(1);
  }
  if(trit == top) {
    std::cerr << "at top can't get parent\n";
    exit(1);
  }
}

// read in an end XML tag
static std::string endXMLtag(std::string & sline){
  std::string::size_type i1 = sline.find("</"); // is it a end key?
  if(i1 == std::string::npos){std::cerr<<"Error! no </ \n";exit(1);}
  std::string::size_type i2 = sline.find(">"); // get tag
  if(i2 == std::string::npos){std::cerr<<"Error! </ but no end > \n";exit(1);}
  std::string sdata = sline.substr(i1+2,i2-i1-2);
  sline.erase(0,i2-i1+1);
  std::string stag;
  std::istringstream ist(sdata); ist >> stag;
  return stag;
}

// parse a single line of input looking for xml tags and end tags comments etc.
std::string & parseline(std::string & sline,
			const std::list<XML>::iterator & top,
			std::list<XML>::iterator & trit,XML & xml1,
			std::stack<std::list<XML>::iterator> & xmls)
{
  std::string::size_type i1 = sline.find("<");
  std::string stag,satt;
  std::string head = HEAD;
  if( i1 == std::string::npos) {
    (*trit).dataadd(sline); 
    sline.erase();
    return sline;
  }
  if(i1>0){
    stag = sline.substr(0,i1); // we found a "<"
    (*trit).dataadd(stag);
    sline = sline.substr(i1);
  }
  if(sline.at(1)=='?') { // is it the header?
    if(sline == head) sline = sline.substr(head.size());
    else  std::cerr << "sline not head? = \"" << sline << "\"\n";
  } else if(sline.at(1)=='/') { // is it an end key?
    stag = endXMLtag(sline);
    errendtag(stag,trit,top);
    trit = xmls.top(); 
    xmls.pop();
    (*trit).dataaddnext("");
    // std::cout << "end tag = " << stag << std::endl;
  } else { // must be a key (ie not comment or end key or head)
    int iend = startXMLtag(sline,stag,satt);
    xml1.tag(stag);  xml1.att(satt);
    // if(trit == top) {std::cerr << "at top \n";exit(1);}
    xmls.push(trit);
    ((*trit).xmllist()).push_back(xml1); // add to list
    trit = ((*trit).xmllist()).end(); // position at the last one!
    --trit; 
    if(iend) { // we also ended the key.. 
      trit = xmls.top();  xmls.pop(); // no data so pop back off
      (*trit).dataaddnext(""); // and go to next string
    }
  }
  return sline;
}

// parses an XML file and returns a list of XML keys parsed
std::list<XML> parseXMLfile(const std::string &filename)
{
  std::list<XML> trx;
  std::list<XML>::iterator top,trit;
  std::stack<std::list<XML>::iterator> xmls;
  XML xml1;

  // std::cout << "Parsing file = \"" << filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if(file.bad()) {
    std::cerr << "Unable to open file \"" << filename << "\"\n";
    std::cerr << "While trying to parse an XML file\n";
    exit(1);
  }
  std::string head = HEAD;
  std::string sline;

  xml1.tag(head); // create head node
  trx.push_back(xml1);  
  trit = top = trx.begin();

  while(getline(file,sline)){
    // std::cout << "in: "<< sline << std::endl;
    sline = skipXMLcomment(file,sline);
    while(sline.size()>0){
      sline = parseline(sline,top,trit,xml1,xmls);
    }
  }

  file.close();

  return trx;
}
