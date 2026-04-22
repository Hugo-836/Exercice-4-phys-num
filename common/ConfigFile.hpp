#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "ConfigFile.h"

/* Constructor of the ConfigFile class
   inputs:
     filename: (str) string containing the file name to read
*/
ConfigFile::ConfigFile(const std::string& fileName){
  // define a ifstream class containing the read file buffer
  std::ifstream file;
  // try to open the file
  // .c_str returns a pointer to an array containing a 
  // null terminated sequence of char terminated '\0'
  file.open(fileName.c_str());
  // i it cannot be opened
  if (!file){
    // send an error message
    std::cerr << "[ConfigFile] Impossible d'ouvrir le fichier " << fileName << std::endl;
  }
  else{
    // init a string for storing input file char
    std::string lineread;
    // until it is possible to read line in files
    while(getline(file, lineread)){
      // process the read input line
      process(lineread);
    }
    // close the input file
    file.close();
  }
}
// destructor of the ConfigFile class 
ConfigFile::~ConfigFile(){}

/* Method for printing data in string format in an output file
   inputs:
     path: (str) string containing the path to the outputfile
*/
void ConfigFile::printOut(const std::string& path) const {
  // open the output file as output stream object
  std::ofstream outputFile(path.c_str());
  // check if the outputfile is open
  if (outputFile.is_open())
  {
    // write ouputdata as strings
    outputFile << toString() << std::endl;
  }
  // close the output file
  outputFile.close();
}

/* Generate a string of the form: string1 = string2 \n */
std::string ConfigFile::toString() const {
  std::string strToReturn; // define a string output

  // Loop on the elements (key,value) within the maps 'configMap'
  for (std::map<std::string,std::string>::const_iterator iter = configMap.begin(); iter != configMap.end(); ++iter) {
      strToReturn.append(iter->first);  // extract key from the map iter element
      strToReturn.append("=");          // add '=' sign
      strToReturn.append(iter->second); // extract value from the map iter element
      strToReturn.append("\n");         // add back to the line
  }
  return strToReturn; // return the string
}

/* The method extract information for a read string add the LHS 
   as key value in configMap map and RHS as value
   inputs:
     lineread: (str) string line read from file
*/
void ConfigFile::process(const std::string& lineread) {
  // find where there is the comment symbol % in the string. Everything at its right
  // then having higher position is considered as comment
  size_t commentPosition=trim(lineread).find('%',0);
  // if there is a key = value element in line (the line is not entirely comment) and the
  // line is not empty (without end of line -> remove tabs and space at both string sides)
  if(commentPosition!=0 && trim(lineread).length()>0){ // End of line is counted as a character on some architectures
    size_t equalPosition = lineread.find('=',1); // find the position of special character '='
    // check if the position of = is not 'valid' (not end of string)
    if(equalPosition==std::string::npos){
      // send an error message
      std::cerr << "[ConfigFile] Ligne sans '=' : \"" << trim(lineread) << "\"" << std::endl;
    }else{
      // identify as key the string substring from the beginning to '=' character (LHS)
      std::string key = trim(lineread.substr(0,equalPosition));
      // identify as value the string substring from after '=' to end (RHS)
      std::string value = trim(lineread.substr(equalPosition+1,lineread.length()));
      // check if there is already a key in the map
      std::map<std::string, std::string>::const_iterator val = configMap.find(key);
      // check if the value is not the end of map
      if (val != configMap.end()) {
        configMap.erase(key); // remove the key
      }
      // add the key with the new value as two strings
      configMap.insert( std::pair<std::string, std::string>(key,value) );
    }
  }
}

/* get procedure: return the value associate to a key with type T
   inputs:  
     key:      (str) map key of the desired value
     initValue (T)(optional) initial value in case the reading 
			     procdedure fails, default: 0/NULL
   outputs:
     out: (T) value with the desired type  
*/
template<typename T> T ConfigFile::get(const std::string& key,const T& initValue) const{
  // extract the iterator of the desired key
  std::map<std::string, std::string>::const_iterator val = configMap.find(key);
  T out(initValue); // output of type T
  // check if the output is valid (not end of map)
  if ( val != configMap.end() ) {
    // initialise input stream as the map value
    std::istringstream iss(val->second);
    iss >> out; // cast value type from input stream
  }else{
    // send an error message
    std::cerr << "[ConfigFile] Le parametre suivant est manquant : " << key << std::endl;
  }
  std::cout << "\t" << key << "=" << out << std::endl; // print the output value
  return out; // return the output
}

/* Retrun a boolean from from key or the current input string state
   if casting fails
   inputs:
     key:      (str) map key of the desired element
     initValue (T)(optional) initial value in case the reading 
			     procdedure fails, default: 0/NULL
   outputs:
     result: (bool) boolean value from input is available of 
                    input string state otherwise
*/
template<> bool ConfigFile::get<bool>(const std::string& key,const bool& initValue) const{
  // extract the value as input string element
  std::istringstream iss(configMap.find(key)->second);
  bool result(initValue); // initialise result as false
  iss >> result; // output the value as boolean
  // check if the conversion has been successful
  if (iss.fail()){
    // set the input string error state with the current fail state
    iss.clear();
    // set result as boolean of the current state
    iss >> std::boolalpha >> result;
  }
  // print result on screen 
  std::cout << "\t" << key << "=" << result << std::endl;
  return result; // return result
}

/* Remove tabs and spaces at the beginning and end of a string
   input:
     str: (str) input string
   output:
     str: (sre) string without beginning / ending tabs and spaces
*/
std::string ConfigFile::trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t"); // find first empty space
    if(first==std::string::npos) // check if the line has only white space
      return ""; // return an empty string
    size_t last = str.find_last_not_of(" \t\r"); // find last white space
    return str.substr(first, (last-first+1)); // return string without white spaces
}
