#ifndef __PARAM_STRING_H__
#define __PARAM_STRING_H__

/** This is a quick-and-dirty set of functions for storing structured
    data in strings.
 */

#include <map>
#include <vector>
#include <string>

using namespace std;

#define PARAM_STRING_FIRST_WORD "Param"

class ParamString {

  map<string,string> params;

 public:
  /** Given a FILE* at the beginning of a parameter section, skip
      ahead to the first byte after the parameter section. Return
      false iff not currently at the beginning of a parameter section */
  static bool skipParameters(FILE *inf);
  
  /** Given a FILE* at the beginning of a parameter section, read
      the parameters into this object. If a syntax error is found,
      print an error to stdout and return false. */
  bool readParameters(FILE *inf);

  /** Write a set of parameters to a file. */
  bool writeParameters(FILE *outf);

  /** Parse a string as an integer. */
  static bool parseInt(const char *s, int &result);

  /** Parse a string as a float. */
  static bool parseFloat(const char *s, float &result);

  /** Parse a string into a list of strings. */
  static bool parseList(const char *s, vector<string> &result);

  /** Parse a string into a set of name/value pairs. */
  static bool parseMap(const char *s, map<string,string> &result);

  static string formatInt(int i);
  static string formatFloat(float f);
  static string formatList(const vector<string> &v);
  static string formatMap(const map<string,string> &m);

  void clear() {params.clear();}
  void set(const char *name, const char *value) {params[name] = value;}

  void setInt(const char *name, int value) {
    params[name] = formatInt(value);
  }

  void setFloat(const char *name, float value) {
    params[name] = formatFloat(value);
  }


  void setList(const char *name, const vector<string> &value) {
    params[name] = formatList(value);
  }

  void setFloatList(const char *name, const vector<float> &floats);

  void setMap(const char *name, const map<string,string> &m) {
    params[name] = formatMap(m);
  }    

  /** Returns true iff a parameter with the given name has been read. */
  bool exists(const char *name) {
    return params.find(name) != params.end();
  }

  /** Get a parameter as a string. Returns false iff the parameter has not
      been read. */
  bool get(const char *name, string &value);

  /** Get a parameter as an integer. Returns false if the parameter has
      not been read, or if it is not parsable as an integer. */
  bool getInt(const char *name, int &value);

  /** Get a parameter as a float. Returns false if the parameter has
      not been read, or if it is not parsable as a float. */
  bool getFloat(const char *name, float &value);

  /** Get a parameter as a list of strings. Returns false if the parameter
      has not been read, or if it is not parsable as a list. */
  bool getList(const char *name, vector<string> &value);

  /** Get a parameter as a list of integers. Returns false if the parameter
      has not been read, or if it is not parsable as a list of integers. */
  bool getIntList(const char *name, vector<int> &value);

  /** Get a parameter as a list of floats. Returns false if the parameter
      has not been read, or if it is not parsable as a list of floats. */
  bool getFloatList(const char *name, vector<float> &value);

  /** Gets a parameter as a set of name/value pairs. */
  bool getMap(const char *name, map<string,string> &value);

 private:
  static const char *getFirstWord() {return "Param";}

  /** Reads a string from an input string, removing escaping.
      This includes lists and name/value maps. For example:
        given: x=y, return "x" and len=1
	given: "x"="y", return "x" and len=3
        given: [123,456]=foo, return "[123,456]" and len=9
  */
  static string readString(const char *&inputStr, const char *endChars);

  enum FormatStringEscapeType {
    FORMAT_STRING_DEFAULT,
    FORMAT_STRING_LHS,
    FORMAT_STRING_RHS
  };

  /** Format a string for serialization. If it contains no specialp
      characters, it is returned unchanged. Otherwise it is quoted,
      and special characters are escaped. */
  static string formatString(const char *s,
			     FormatStringEscapeType esc=FORMAT_STRING_DEFAULT);

  /** Returns a pointer to the first non-whitespace character after
      the given pointer. */
  static const char *skipWhitespace(const char *s);

  /** Read the first line in a parameter section. Returns the length
      int bytes of the rest of the section, or a negative number
      on error. */
  static int readFirstLine(FILE *inf);

  /** Read one line of data from the input file. Return false on EOF. */
  static bool readLine(FILE *inf, string &line);

};


#endif // __PARAM_STRING_H__
