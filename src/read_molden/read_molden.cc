#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#define Ang2Bohr 1.889725989

using namespace std;

int main(int argc, char *argv[]) {
   if (argc != 2){
      cout << "No filename present!" << endl;
      return -1;
   }

   string filename(argv[1]);

   ifstream ifs;
   ifs.open(filename.c_str());

   /************************************************************
   *  Set up variables that will contain the organized info    *
   ************************************************************/
   int num_basis = 0, num_atoms = 0;

   /* Atom positions */
   vector< vector<double> > positions;
   /* Atom names */
   vector< string > names;
   /* Vector of atomic basis information vectors */
   vector<vector<tuple<string, vector<double>, vector<vector<double> > > > > basis_info;
   /* atom order in the gto section */
   vector<int> gto_order;
   /* Matrix of coefficients... not yet fully implemented */
   vector< vector<double> > mo_coefficients;

   /************************************************************
   *  Set up "global" regular expressions                      *
   ************************************************************/
   boost::regex gto_re("\\[GTO\\]");
   boost::regex atoms_re("\\[Atoms\\]");
   boost::regex mo_re("\\[MO\\]");
   boost::regex other_re("\\[\\w+\\]");

   boost::cmatch matches;

   /************************************************************
   *  Booleans to check and make sure each important section   *
   *  was found                                                *
   ************************************************************/
   bool found_atoms = false;
   bool found_gto = false;
   bool found_mo = false;

   double scale;

   string line; // Contains the current line of the file
   getline(ifs, line);

   while (!ifs.eof()){
      if (boost::regex_search(line,atoms_re)) {
         boost::regex ang_re("Angs");
         boost::regex atoms_line("(\\w{1,2})\\s+\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");

         scale = boost::regex_search(line, ang_re) ? Ang2Bohr : 1.0;

         getline(ifs, line);
         while(!boost::regex_search(line, other_re)){
            if (ifs.eof()) { break; }

            if(boost::regex_search(line.c_str(), matches, atoms_line)) {
               ++num_atoms;

               string nm(matches[1].first, matches[1].second);
               names.push_back(nm);
               
               string x_str(matches[2].first, matches[2].second);
               string y_str(matches[3].first, matches[3].second);
               string z_str(matches[4].first, matches[4].second);

               vector<double> pos;

               pos.push_back(boost::lexical_cast<double>(x_str));
               pos.push_back(boost::lexical_cast<double>(y_str));
               pos.push_back(boost::lexical_cast<double>(z_str));
               
               positions.push_back(pos);

               getline(ifs,line);
            }

            else { getline(ifs,line); }
         }

         found_atoms = true;
      }
      else if (boost::regex_search(line,gto_re)){
         getline(ifs, line);

         boost::regex atom_line("(\\d+)\\s+\\d?");
         boost::regex shell_line("([spd])\\s+(\\d+)\\s+\\S*");
         boost::regex exp_line("(\\S+)\\s+(\\S+)");
         boost::regex Dd("[Dd]");

         while(!boost::regex_search(line,other_re)){
            if (ifs.eof()) { break; }

            /* This line should be a new atom */
            if(!boost::regex_search(line.c_str(), matches, atom_line)) {
               getline(ifs,line); continue;
            }

            vector<double>::iterator last_set_start;
            vector<double>::iterator last_set_end;

            int offset = 0;

            vector<double> exponents;
            vector< vector< double > > coefficients;
            vector<tuple<string,vector<double>,vector<vector<double> > > > atom_basis_info;

            string atom_no(matches[1].first, matches[1].second);
            gto_order.push_back(boost::lexical_cast<int>(atom_no.c_str()));

            string ang_l;

            getline(ifs, line);

            while(boost::regex_search(line.c_str(), matches, shell_line)) {
               /* Now it should be a new angular shell */
               string tmp_ang_l(matches[1].first, matches[1].second);
               if (tmp_ang_l != ang_l) { // new l value, new shell
                  // package old shell data (tuple including <ang_l, exponents, coefficients>)
                  atom_basis_info.push_back(make_tuple(ang_l, exponents, coefficients));
                  // initialize new shell data
                  exponents.clear();
                  coefficients.clear();
                  ang_l = tmp_ang_l;
                  // reset offset
                  offset = 0;
               }
               // TODO: figure out how to read and include [5D] keyword
               if(ang_l == "s")        { num_basis += 1; }
               else if(ang_l == "p")   { num_basis += 3; }
               else if(ang_l == "d")   { num_basis += 6; }

               string num_exp_string(matches[2].first, matches[2].second);
               int num_exponents = boost::lexical_cast<int>(num_exp_string);

               vector<double> evec;
               vector<double> cvec;

               //cout << ang_l << "   " << num_exponents << endl;

               for (int i = 0; i < num_exponents; ++i) {
                  getline(ifs,line);

                  boost::regex_search(line.c_str(), matches, exp_line);

                  string exp_string(matches[1].first, matches[1].second);
                  string coeff_string(matches[2].first, matches[2].second);

                  exp_string = boost::regex_replace(exp_string, Dd, "E");
                  coeff_string = boost::regex_replace(coeff_string, Dd, "E");

                  double exponent = boost::lexical_cast<double>(exp_string);
                  double coeff = boost::lexical_cast<double>(coeff_string);

                  evec.push_back(exponent);
                  cvec.push_back(coeff);

                  //cout << exponent << "   " << coeff << endl;
               }

               /************************************************************
               *  Right now I have evec and cvec which have the newly      *
               *  exponents and coefficients                               *
               ************************************************************/
               if( !exponents.empty() ) {
                  bool same_set;
                  if ( distance(last_set_start, last_set_end) == num_exponents ) {
                     if (equal (last_set_start, last_set_end, evec.begin() ) ) { same_set = true; }
                     else { same_set = false; }
                  }
                  else { same_set = false; }
            
                  if(!same_set) {
                     offset += distance(last_set_start, last_set_end);

                     // insert new exponents
                     exponents.insert(exponents.end(), evec.begin(), evec.end());

                     // update iterators
                     last_set_end   = exponents.end();
                     last_set_start = last_set_end - num_exponents;
                  }
               }
               else {
                  // insert new exponents
                  exponents.insert(exponents.end(), evec.begin(), evec.end());

                  last_set_start = exponents.begin();
                  last_set_end = exponents.end();
               }

               // add coefficients for the new contracted ao
               vector<double> tmp_cvec(offset, 0.0);
               tmp_cvec.insert(tmp_cvec.end(), cvec.begin(), cvec.end());
               coefficients.push_back(tmp_cvec);

               getline(ifs, line); 
            }
            // add current shell to atom_basis_info
            atom_basis_info.push_back(make_tuple(ang_l, exponents, coefficients));
            //atom_basis_info to basis_info
            basis_info.push_back(atom_basis_info);

            //cout << endl;
         }

         found_gto = true;
      }
      else if (boost::regex_search(line,mo_re)) {
         break; // this section needs work
         if(found_gto) {
            //cout << "num_basis = " << num_basis << endl;
         }
         else {
            cout << "GTO section hasn't been found yet!" << endl; break;
         }

         boost::regex ene_re("Ene=\\s+(\\S+)");
         boost::regex spin_re("Spin=\\s+(\\w+)");
         boost::regex occup_re("Occup=\\s+(\\S+)");
         boost::regex coeff_re("\\d+\\s+(\\S+)");

         while(!boost::regex_search(line,other_re)) {
            if (ifs.eof()) { break; }

            vector<double> movec;

            /* MO Energy */
            getline(ifs,line);
            if(!boost::regex_search(line.c_str(), matches, ene_re)) {
               cout << "Whoops, there's a problem in with the MO section" << endl;
            }

            /* Spin */
            getline(ifs,line);
            if(!boost::regex_search(line.c_str(), matches, spin_re)) {
               cout << "Whoops, there's a problem in with the MO section" << endl;
            }
            
            /* Occupation */
            getline(ifs,line);
            if(!boost::regex_search(line.c_str(), matches, occup_re)) {
               cout << "Whoops, there's a problem in with the MO section" << endl;
            }

            /* Read in MO coefficients */
            for (int i = 0; i < num_basis; ++i) {
               getline(ifs,line);
               boost::regex_search(line.c_str(), matches, coeff_re);
               string mo_string(matches[1].first, matches[1].second);
               double coeff = boost::lexical_cast<double>(mo_string);
               cout << i+1 << "   " << coeff << endl;
               movec.push_back(coeff);
            }
            mo_coefficients.push_back(movec);
         }
      }
      else {
         getline(ifs,line);
      }
   }

   ifs.close();

   vector< string >::iterator niter = names.begin();
   for (vector< vector< double > >::iterator piter = positions.begin(); piter != positions.end(); ++piter, ++niter) {
      cout << *niter << "   " << (*piter)[0] << "   " << (*piter)[1] << "   " << (*piter)[2] << endl;
   }
   cout << endl << endl;

   /* At this point, I hypothetically have a vector of vectors of all the info I need for the shells. Now to print. */
   /* This prints everything in atom order, i.e., gto order doesn't have to match atom order */
   for(int i = 0; i < num_atoms; ++i){
      auto bviter = basis_info.begin();
      for(vector<int>::iterator oiter = gto_order.begin(); oiter != gto_order.end(); ++oiter) {
         if( (i+1) == *oiter ) { break; }
         else { ++bviter; }
      }

      cout << "new atom line";
      for(auto aiter = bviter->begin(); aiter != bviter->end(); ++aiter) {
         const string current_name = get<0>(*aiter);
         cout << current_name << endl;

         const vector<double> current_exp = get<1>(*aiter);
         const vector<vector<double> > current_coeff = get<2>(*aiter);
         
         int ii = 0;

         for(vector<double>::const_iterator eiter = current_exp.begin(); eiter != current_exp.end(); ++eiter) {
            cout << *eiter;

            for (vector<vector<double> >::const_iterator citer = current_coeff.begin(); citer != current_coeff.end(); ++citer) {
               if( ii < citer->size() ) {
                  cout << "   " << (*citer)[ii];
               }
               else {
                  cout << "   " << "--";
               }
            }
            ++ii;
            cout << endl;
         }
         cout << endl;
      }
   }

   /* TODO: Add in print statement to check mo_coefficients */

   return 0;
}
