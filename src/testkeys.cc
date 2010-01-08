
#include <stdlib.h>
#include <algorithm>
#include <list>
#include "simkeys.h++"

int main()
{
  SIMKEY tkey,nkey;
  tkey.key("happy","days");
  nkey = tkey;
  if(nkey==tkey) cout << "keys are equil" << endl;
  else cout << "keys aren't equil" << endl;
  nkey.key("new","fellow");
  if(nkey==tkey) cout << "keys are equil" << endl;
  else cout << "keys aren't equil" << endl;
  list<SIMKEY> keys;
  list<SIMKEY>::iterator iter;
  keys.push_back(tkey);
  keys.push_back(nkey);

  tkey.key("new");
  iter=find(keys.begin(),keys.end(), "new");
  cout << *iter <<endl;
  return 0;
}
