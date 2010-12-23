/*
 * main.cpp
 *
 *  Created on: Oct 20, 2010
 *      Author: zhangz
 */

#include <flann/flann.hpp>
#include <flann/io/hdf5.h>
#include <dirent.h>
#include <stdio.h>

#include <time.h>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;

/*
 * input
 * 	dir:	full path to directory containing sift files
 * output
 * 	files:	vector of sift filenames in directory
 * returns
 * 	0 if successful, -1 if error
 */
int getSIFTFilenames(string dir, vector<string> &files) {
	DIR *dp;
	struct dirent *dirp;
	//cout<<dir<<endl;
	if ((dp = opendir(dir.c_str())) == NULL) {
		cout << "Error opening " << dir << endl;
		return -1;
	}
	while ((dirp = readdir(dp)) != NULL) {
		//if (dirp->d_type == DT_REG) {
		string fname = string(dirp->d_name);
		if (fname.find("sift.txt") != string::npos)
			files.push_back(fname);
		//}
	}
	closedir(dp);
	return 0;
}
/*
 * input
 * 	files:	vector of sift filenames
 * 	dir:	directory files are located in
 * returns
 * 	total number of features in files, -1 if error
 */
int getNumFeaturesInFile(string filename) {
	ifstream file;
	file.open(filename.c_str(), ios::in);
	if (file.is_open()) {
		string line;
		getline(file, line);
		file.close();
		return atoi(line.substr(0, line.find(' ')).c_str());
	} else {
		cout << "error opening file: " << filename << endl;
		return -1;
	}
}
/*
 * input
 * 	files:	vector of sift filenames
 * 	dir:	directory files are located in
 * returns
 * 	total number of features in files, -1 if error
 */
int getNumFeaturesInDir(string dir, vector<string> &files) {
	ifstream file;
	unsigned int numfeats = 0;
	for (unsigned int i = 0; i < files.size(); i++) {
		string fname = files[i];
		string path = dir + '/' + fname;
		file.open(path.c_str(), ios::in);
		if (file.is_open()) {
			string line;
			getline(file, line);
			numfeats += atoi(line.substr(0, line.find(' ')).c_str());
		} else {
			cout << "error opening file: " << path << endl;
			return -1;
		}
		file.close();
	}
	return numfeats;
}

/*
 * input
 * 	file:	path of file to write to
 * 	vec:	vector of strings to write out to file
 * returns
 * 	0 if successful, -1 if fail
 */
int saveVectorToFile(string file, vector<string> vec) {
	std::ofstream outfile(file.c_str());
	if (!outfile.is_open())
		return -1;
	for (unsigned int i = 0; i < vec.size(); i++) {
		outfile.write(vec[i].c_str(), vec[i].length());
		outfile.put('\n');
	}
	outfile.close();
	return 0;
}
/*
 * input
 * 	file:	path of file to read from
 * 	vec:	vector of strings to read file into
 * returns
 * 	0 if successful, -1 if fail
 */
int readVectorFromFile(string file, vector<string> &vec) {
	std::ifstream infile(file.c_str());
	if (!infile.is_open())
		return -1;
	std::string line;
	while (getline(infile, line)) {
		vec.push_back(line);
	}
	infile.close();
	return 0;
}

int readInSiftFile(string filepath, flann::Matrix<int> &features) {
	int row = 0;
	std::ifstream myfile(filepath.c_str(), ios::in);
	if (myfile.is_open()) {
		string line;
		getline(myfile, line); //ignore file header
		while (myfile.good()) {
			//read in a feature
			getline(myfile, line); //ignore feature header
			if (myfile.eof())
				break;
			int col = 0;
			for (int x = 0; x < 6; x++) {
				getline(myfile, line, ' ');
				for (int y = 0; y < 19; y++) {
					getline(myfile, line, ' ');
					features[row][col] = atoi(line.c_str());
					col++;
				}
				getline(myfile, line, '\n');
				features[row][col] = atoi(line.c_str());
				col++;
			}
			getline(myfile, line, ' ');
			for (int y = 0; y < 7; y++) {
				getline(myfile, line, ' ');
				features[row][col] = atoi(line.c_str());
				col++;
			}
			getline(myfile, line, '\n');
			features[row][col] = atoi(line.c_str());
			col++;
			row++;
		}
	} else {
		cout << "error opening file: " << filepath << endl;
	}
	//cout<<"finished reading file: "<<fname<<endl;
	myfile.close();
	return 0;
}
int readInSiftFiles(string dir, vector<string> &files,
		flann::Matrix<int> &dataset, vector<string> &feat_to_img) {
	int row = 0;
	for (unsigned int i = 0; i < files.size(); i++) {
		string fname = files[i];
		string path = dir + '/' + fname;
		std::ifstream myfile(path.c_str(), ios::in);
		if (myfile.is_open()) {
			string line;
			getline(myfile, line); //ignore file header
			while (myfile.good()) {
				//read in a feature
				getline(myfile, line); //ignore feature header
				if (myfile.eof())
					break;
				int col = 0;
				for (int x = 0; x < 6; x++) {
					getline(myfile, line, ' ');
					for (int y = 0; y < 19; y++) {
						getline(myfile, line, ' ');
						dataset[row][col] = atoi(line.c_str());
						col++;
					}
					getline(myfile, line, '\n');
					dataset[row][col] = atoi(line.c_str());
					col++;
				}
				getline(myfile, line, ' ');
				for (int y = 0; y < 7; y++) {
					getline(myfile, line, ' ');
					dataset[row][col] = atoi(line.c_str());
					col++;
				}
				getline(myfile, line, '\n');
				dataset[row][col] = atoi(line.c_str());
				col++;
				feat_to_img.push_back(files[i]);//fname;
				row++;
			}
		} else {
			cout << "error opening file: " << path << endl;
		}
		myfile.close();
	}
	return 0;
}
bool fexists(string filename) {
	ifstream ifile(filename.c_str());
	return ifile;
}

void printer(multimap<float, string> pN) {
	cout << "Map size = " << pN.size() << endl;
	multimap<float, string>::iterator it = pN.begin();
	int count = 0;
	while (it != pN.end()) {
		count += it->first;
		cout << "Key = " << it->first << "\t\t Value = " << it->second << endl;
		it++;
	}
	cout << count << endl;
}

int alg(string celldir, string cell, string qdir, string qfile,
		string *outfile = NULL, string indextype = "kdtree1", bool output =
				false) {
	string qpath = qdir + qfile;
	string cellpath = celldir + cell;

	flann::Matrix<int> dataset;
	vector<string> feat_to_img = vector<string> ();

	time_t start_time = time(NULL);
	time_t last_op_time = start_time;

	cout << "query file: " << qfile << endl;
	cout << "cell: " << cell << endl;

	//read in the cell dataset
	if (fexists(cellpath + "dataset.hdf5") && fexists(cellpath
			+ "feat_to_img.dat")) {
		cout << "loading db from file" << endl;
		flann::load_from_file(dataset, cellpath + "dataset.hdf5", cell
				+ "dataset");
		readVectorFromFile(cellpath + "feat_to_img.dat", feat_to_img);
		cout << "db loaded" << endl;
	} else {

		cout << "finding sift files" << endl;
		//get sift files in dir
		vector<string> files = vector<string> ();
		if (getSIFTFilenames(cellpath, files) != 0) {
			return 1;
		}

		//gets number of features in dir
		dataset.cols = 128; //TODO: put constants in one place
		if ((dataset.rows = getNumFeaturesInDir(cellpath, files)) <= 0) {
			return 1;
		}
		dataset.data = new int[dataset.rows * dataset.cols];
		//flann::Matrix<string> feat_to_img(new string[dataset.rows], dataset.rows, 1);

		//		cout << "reading in sift files" << endl;
		readInSiftFiles(cellpath, files, dataset, feat_to_img);

		//cout << "saving data" << endl;

		//flann::save_to_file(dataset, cellpath + "dataset.hdf5", "dataset");

		//saveVectorToFile(cellpath + "feat_to_img.dat", feat_to_img);
		//cout << "data saved" << endl;
	}
	cout << "took " << time(NULL) - last_op_time << " seconds" << endl;
	cout << "features in db: " << dataset.rows << endl;
	last_op_time = time(NULL);

	//read in the query dataset
	flann::Matrix<int> query;
	if (fexists(qpath + "query.hdf5")) {
		cout << "loading query from file" << endl;
		flann::load_from_file(query, qpath + "query.hdf5", qfile + "query");
		cout << "query loaded" << endl;
	} else {
		cout << "reading in query" << endl;
		query.cols = 128; //TODO: put constants in one place
		query.rows = getNumFeaturesInFile(qpath);
		query.data = new int[query.rows * query.cols];
		readInSiftFile(qpath, query);
		cout << "saving query" << endl;
		flann::save_to_file(query, qpath + "query.hdf5", qfile + "query");
	}

	cout << "took " << time(NULL) - last_op_time << " seconds" << endl;
	cout << "features in query: " << query.rows << endl;
	last_op_time = time(NULL);

	string indexsuffix = indextype + ".index";

	flann::Index<flann::L2<int> > *index;
	if (fexists(cellpath + indexsuffix)) {
		cout << "loading index" << endl;
		index = new flann::Index<flann::L2<int> >(dataset, flann::SavedIndexParams(cellpath
				+ indexsuffix));
	} else {
		cout << "creating index" << endl;
		if (indextype == "kdtree4")
			// construct an randomized kd-tree index using 4 kd-trees
			index = new flann::Index<flann::L2<int> >(dataset, flann::KDTreeIndexParams(4));
		else if (indextype == "kdtree1")
			index = new flann::Index<flann::L2<int> >(dataset, flann::KDTreeIndexParams(1));
		else if (indextype == "kmeans64,22")
			index = new flann::Index<flann::L2<int> >(dataset, flann::KMeansIndexParams(64,
					22));
		else {
			cout << "invalid index specified" << endl;
			return -1;
		}
		index->buildIndex();
		index->save(cellpath + indexsuffix);
	}

	cout << "took " << time(NULL) - last_op_time << " seconds" << endl;
	last_op_time = time(NULL);
	// do a knn search, using 128 checks
	cout << "querying index" << endl;
	int nn = 1;
	flann::Matrix<int> indices(new int[query.rows * nn], query.rows, nn);
	flann::Matrix<float> dists(new float[query.rows * nn], query.rows, nn);
	index->knnSearch(query, indices, dists, nn, flann::SearchParams(2048));
	//votes
	map<string, float> counts;
	map<string, int> sums;

	//	ofstream distfile;
	//	distfile.open("distances.lst");
	//	float a = 50001210;
	for (unsigned int i = 0; i < dists.rows; i++) {
		if (dists[i][0] < 70000) {
			counts[feat_to_img[indices[i][0]]] += 1;
			//counts[feat_to_img[indices[i][0]]] += a / (a + (dists[i][0]));
			//counts[feat_to_img[indices[i][0]]] += a - dists[i][0];
		}
		//distfile<<dists[i][0]<<endl;
	}
	//	distfile.close();


	multimap<float, string> pN;

	cout << "took " << time(NULL) - last_op_time << " seconds" << endl;
	last_op_time = time(NULL);

	cout << "total time of " << time(NULL) - start_time << " seconds" << endl;

	map<string, float>::iterator countitr;
	for (countitr = counts.begin(); countitr != counts.end(); countitr++)
		pN.insert(pair<float, string> (countitr->second, countitr->first));

	if (output)
		printer(pN);

	if (outfile != NULL) {
		ofstream myfile;
		myfile.open(outfile->c_str());
		multimap<float, string>::reverse_iterator it = pN.rbegin();
		while (it != pN.rend()) {
			myfile << it->first << "\t" << it->second << endl;
			it++;
		}
		myfile.close();
	}

	delete index;
	dataset.free();
	query.free();
	indices.free();
	dists.free();

	return 0;

}

//bool comp()

int main(int argc, char** argv) {

	if (argc == 5) {
		string celldir = string(argv[1]); //path to directory containing cells
		string cell = string(argv[2]); //name of cell (directory)
		string qdir = string(argv[3]); //path to directory containing query
		string qfile = string(argv[4]); //name of query (sift file)
		alg(celldir, cell, qdir, qfile);
	} else if (argc == 6) {
		string celldir = string(argv[1]); //path to directory containing cells
		string cell = string(argv[2]); //name of cell (directory)
		string qdir = string(argv[3]); //path to directory containing query
		string qfile = string(argv[4]); //name of query (sift file)
		string outfile = string(argv[5]); //path of outputfile
		alg(celldir, cell, qdir, qfile, &outfile);
	} else if (argc == 7) {
		string celldir = string(argv[1]); //path to directory containing cells
		string cell = string(argv[2]); //name of cell (directory)
		string qdir = string(argv[3]); //path to directory containing query
		string qfile = string(argv[4]); //name of query (sift file)
		string outfile = string(argv[5]); //path of outputfile
		string index = string(argv[6]); //searchindex
		alg(celldir, cell, qdir, qfile, &outfile, index);
	}
}

