/*
 * objReader.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: zxi
 */
#include "objReader.h"

#include "util/Path.h"


//convert a string to a list of tokens
inline list<string> tokenize(char * tmp, const char * ignore)
{
	list<string> tokens;
	char * tok = strtok(tmp, ignore);
	while (tok != NULL) {
		tokens.push_back(tok);
		tok = strtok(NULL, ignore);
	}
	return tokens;
}

inline void getAllLabels(istream& in, list<list<string> >& tokens)
{
	const int size = 1024;
	char * tmp = new char[size];

	while (!in.eof()) {
		in.getline(tmp, size);
		//check for termination
		if (string(tmp).find("}") != string::npos) //found "}"
			break;

		list<string> tok = tokenize(tmp, " \t[]()<>,={");
		if (tok.empty())
			continue;
		string label = tok.front();
		if (label[0] == '#')
			continue; //comment

		tokens.push_back(tok);
	} //end while

	delete[] tmp;
}


istream& operator>>(istream& in, V& v) {
  in >> v.x >> v.y >> v.z;
  return in;
}

ostream& operator<<(ostream& out, const V& v) {
  out << v.x << " " << v.y << " " << v.z;
  return out;
}

bool objReader::Read(istream& in)
{
  list<list<string> > tokens;
  getAllLabels(in, tokens);

  //read pts
  for( list<string> & token : tokens)
  {
    string label=token.front();
    token.pop_front();
    if (label == "v"){
      vector<double> values;
      for(string& tmp : token) values.push_back(atof(tmp.c_str()));
      Vpt pt;
      pt.x=values[0];
      pt.y=values[1];
      pt.z=values[2];
      m_data.pts.push_back(pt);
      if(values.size()>3)
      {
        Vpt pt;
        pt.x=values[3];
        pt.y=values[4];
        pt.z=values[5];
        m_data.colors.push_back(pt);
      }
    }
    else if (label == "vn"){
      vector<double> values;
      for(string& tmp : token) values.push_back(atof(tmp.c_str()));
      V pt;
      pt.x=values[0];
      pt.y=values[1];
      pt.z=values[2];
      m_data.normals.push_back(pt);
    }
    else if (label == "vt"){
      vector<double> values;
      for(string& tmp : token) values.push_back(atof(tmp.c_str()));
      V pt;
      pt.x=values[0];
      pt.y=values[1];
      m_data.textures.push_back(pt);
    }
    else if(label == "usemtl")
    {
        //TODO
    }
    else if (label == "mtllib")
    {
        this->m_mtl_path=token.front();
        //cout << "mtllib " << mtl_path << std::endl;
        //this->m_mtl_path = mtl_path;
        // assuming obj and mtl are in the same dir.
        string full_path = masc::path::DirPath(this->m_filename) + "/" + this->m_mtl_path;
        this->m_mtl_reader.Read(full_path);
    }
    else if (label == "f")
    {
      polygon poly;

      for(string& tmp : token) //each token value
      {
        int pos1 = (int)tmp.find('/');
        int pos2 = (int)tmp.rfind('/');

        int id_v = std::atoi(tmp.substr(0, pos1).c_str()) - 1;
        int id_t = std::atoi(tmp.substr(pos1 + 1, pos2).c_str()) - 1;
        int id_n = std::atoi(tmp.substr(pos2 + 1).c_str()) - 1;

        poly.pts.push_back(id_v);
        poly.normals.push_back(id_n);
        poly.textures.push_back(id_t);

        // string field1 = tmp.substr(0, pos1);
        // string field2 = tmp.substr(pos1 + 1, pos2 - pos1 - 1);
        // string field3 = tmp.substr(pos2 + 1);
        //
        // if (pos1 < 0 || pos2 < 0) //has no "/"
        // {
        //   field2.clear();
        //   field3.clear();
        // }
        // else if (pos1 == pos2 && pos1 >= 0) //has only on "/"
        // {
        //   field3.clear();
        // }
        //
        // int id_v = atoi(field1.c_str()) - 1;
        // poly.pts.push_back(id_v);
        //
        //
        // if (field2.empty() == false)
        // {
        //   int id_t = atoi(field2.c_str()) - 1;
        //   if (id_t >= 0) poly.textures.push_back(id_t);
        // }
        //
        // if (field3.empty() == false)
        // {
        //   int id_n = atoi(field3.c_str()) - 1;
        //   if (id_n >= 0) poly.normals.push_back(id_n);
        // }
      }//end parsing face

      if (!poly.pts.empty()) m_data.polys.push_back(poly);
    }
    else{}
  }

  m_data.compute_v_normal();

  return true;
}
