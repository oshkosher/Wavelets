#include <cstdio>
#include <cstring>
#include <cctype>
#include "param_string.h"

bool ParamString::skipParameters(FILE *inf) {
  int len = readFirstLine(inf);
  if (len < 0) return false;

  return 0 == fseek(inf, len, SEEK_CUR);
}

bool ParamString::readParameters(FILE *inf) {

  readFirstLine(inf);

  params.clear();
  string line, name, value;
  
  while (true) {
    if (!readLine(inf, line)) return false;
    if (line.length() == 0) break;
    
    const char *str = line.c_str();
    // extract the name--everything until the equals sign
    name = readString(str, "=");
    // skip this line if there is no equals or no name
    if (name.length() == 0 || *str != '=') continue;
    str++;
    // extract the value--everything until the end of the line
    value = readString(str, "");
    params[name] = value;
  }
  return true;
}


bool ParamString::writeParameters(FILE *outf) {
  string content;
  map<string,string>::const_iterator iter;
  for (iter = params.cbegin(); iter != params.cend(); iter++) {
    content.append(formatString(iter->first.c_str(), FORMAT_STRING_LHS));
    content.push_back('=');
    content.append(formatString(iter->second.c_str(), FORMAT_STRING_RHS));
    content.push_back('\n');
  }
  content.push_back('\n');
  size_t size = content.length();
  if (fprintf(outf, PARAM_STRING_FIRST_WORD " %llu\n", (long long unsigned) size) < 0) return false;
  if (fputs(content.c_str(), outf) == EOF) return false;
  return true;
}



bool ParamString::parseInt(const char *s, int &result) {
  return 1 == sscanf(s, "%d", &result);
}


bool ParamString::parseFloat(const char *s, float &result) {
  return 1 == sscanf(s, "%f", &result);
}


bool ParamString::parseList(const char *s, vector<string> &result) {
  result.resize(0);

  s = skipWhitespace(s);
  if (*s != '[') return false;
  s++;
  
  while(true) {
    string word = readString(s, ",]");
    if (*s == '\0') return false;
    result.push_back(word);
    if (*s == ']') break;
    s++;
  }

  return true;
}    


bool ParamString::parseMap(const char *s, map<string,string> &result) {
  result.clear();

  s = skipWhitespace(s);
  if (*s != '{') return false;
  s++;
  string name, value;
  
  while(true) {
    name = readString(s, "=");
    if (*s == '\0') return false;
    s++;
    value = readString(s, ",}");
    if (*s == '\0') return false;

    result[name] = value;

    if (*s == '}') break;
    s++;
  }

  return true;
}    
  

string ParamString::formatInt(int i) {
  char buf[20];
#ifdef _WIN32
  sprintf_s(buf, (sizeof buf) - 1, "%d", i);
#else
  snprintf(buf, sizeof buf, "%d", i);
#endif
  return string(buf);
}

string ParamString::formatFloat(float f) {
  char buf[20];
#ifdef _WIN32
  sprintf_s(buf, (sizeof buf) - 1, "%.8g", f);
#else
  snprintf(buf, sizeof buf, "%.8g", f);
#endif

  return string(buf);
}

string ParamString::formatList(const vector<string> &v) {
  string result;
  result.push_back('[');
  for (size_t i=0; i < v.size(); i++) {
    if (i > 0) result.push_back(',');
    result.append(formatString(v[i].c_str()));
  }
  result.push_back(']');
  return result;
}

string ParamString::formatMap(const map<string,string> &m) {
  string result;
  map<string,string>::const_iterator iter;
  bool first = true;
  result.push_back('{');

  for (iter = m.cbegin(); iter != m.cend(); iter++) {
    if (!first) result.push_back(',');
    first = false;
    result.append(formatString(iter->first.c_str()));
    result.push_back('=');
    result.append(formatString(iter->second.c_str()));
  }

  result.push_back('}');
  return result;
}


void ParamString::setFloatList(const char *name, const vector<float> &floats) {
  vector<string> strings(floats.size());
  for (size_t i=0; i < floats.size(); i++)
    strings[i] = formatFloat(floats[i]);
  setList(name, strings);
}


bool ParamString::get(const char *name, string &value) {
  map<string,string>::iterator p = params.find(name);
  if (p == params.end()) return false;
  value = p->second;
  return true;
}


bool ParamString::getInt(const char *name, int &value) {
  string s;
  if (!get(name, s)) return false;
  return parseInt(s.c_str(), value);
}


bool ParamString::getFloat(const char *name, float &value) {
  string s;
  if (!get(name, s)) return false;
  return parseFloat(s.c_str(), value);
}


bool ParamString::getList(const char *name, vector<string> &value) {
  string s;
  if (!get(name, s)) return false;
  return parseList(s.c_str(), value);
}


bool ParamString::getFloatList(const char *name, vector<float> &values) {
  vector<string> strings;
  if (!getList(name, strings)) return false;
  values.resize(strings.size());
  for (size_t i=0; i < strings.size(); i++) {
    float f;
    if (!parseFloat(strings[i].c_str(), f)) return false;
    values[i] = f;
  }
  return true;
}


bool ParamString::getMap(const char *name, map<string,string> &values) {
  string s;
  if (!get(name, s)) return false;
  return parseMap(s.c_str(), values);
}


/** If the string starts with a doublequote, read until the matching
    doublequote. Otherwise, read until the end of the string, or until
    one of the characters in endChars is reached.
    Advances 'inputStr' by the number of characters consumed.
*/
string ParamString::readString(const char *&inputStr, const char *endChars) {
  string result;
  const char *p = skipWhitespace(inputStr);

  // read until matching doublequote
  if (*p == '"') {
    p++;
    while (*p && *p != '"') {
      if (*p == '\\') {
	switch (p[1]) {
	case '\\': result.push_back('\\'); p++; break;
	case '"': result.push_back('='); p++; break;
	case 'n': result.push_back('\n'); p++; break;
	case 'r': result.push_back('\r'); p++; break;
	case 't': result.push_back('\t'); p++; break;
	default:
	  // not a recognized escape, leave the backslash
	  result.push_back(*p);
	}
      } else {
	result.push_back(*p);
      }
      p++;
    }
    if (*p == '"') p++;
    inputStr = p;
  }

  // read until an endChar
  else {
    const char *e = strpbrk(p, endChars);
    if (e == NULL)
      e = p + strlen(p);
    result.append(p, e);

    // trim trailing whitespace
    while (result.length() > 0 && isspace(result[result.length()-1]))
      result.resize(result.length()-1);

    inputStr = e;
  }

  return result;
}


string ParamString::formatString(const char *s, FormatStringEscapeType esc) {

  // enclose in doublequotes if it contains unprintables or
  // param-file punctuation characters

  const char *endChars;
  if (esc == FORMAT_STRING_LHS) endChars = "\"\n\t =";
  else if (esc == FORMAT_STRING_RHS) endChars = "\"\n\t ";
  else endChars = "\"\n\t ,[]{}=";

  if (!strpbrk(s, endChars)) return string(s);

  // escape doublequotes, newlines, and tab characters
  // for the rest, it's sufficient to enclose them in doublequotes
  string result;
  const char *p = s;
  result.push_back('"');
  while (*p) {
    if (*p == '"') {
      result.append("\\\"");
    } else if (*p == '\n') {
      result.append("\\n");
    } else if (*p == '\t') {
      result.append("\\t");
    } else {
      result.push_back(*p);
    }
    p++;
  }
  result.push_back('"');

  return result;
}

const char *ParamString::skipWhitespace(const char *s) {
  while (isspace(*s)) s++;
  return s;
}


int ParamString::readFirstLine(FILE *inf) {
  string line;
  int len;

  if (!readLine(inf, line)) return -1;
  if (1 != sscanf(line.c_str(), PARAM_STRING_FIRST_WORD " %d", &len))
    return -2;

  if (len < 0) return -3;

  return len;
}

bool ParamString::readLine(FILE *inf, string &line) {
  line.clear();
  int c;
  while (true) {
    c = getc(inf);
    if (c == EOF || c == '\n') break;
    line.push_back((char)c);
  }
  return c != EOF || !line.empty();
}

