//
// Created by Zakhar on 16/03/2017.
//

#include <iostream>
#include <dirent.h>

#include "MatrixReadWrite.h"
#include "../Algebra/MissingBlock.hpp"

using namespace std;

namespace MathIO
{

//
// Reading
//

// global vars
ifstream ___global_file;
string ___global_buffer;

double *getMatrixByName(string name, uint64_t n, uint64_t m)
{
    ___global_file = ifstream(name);
    
    if (!___global_file.is_open())
    {
        std::cout << "Can't open file " << name << std::endl;
        return nullptr;
    }
    
    double *data = (double *)malloc(n * m * sizeof(*data));
    
    uint64_t i = 0;
    
    while (!___global_file.eof() && i < n)
    {
        getline(___global_file, ___global_buffer);
        
        size_t pos = 0;
        
        string temp;
        for (uint64_t j = 0; j < m; ++j)
        {
            // TBA
            size_t newpos;
            newpos = ___global_buffer.find(' ', pos + 1);
            
            newpos = newpos == string::npos ? ___global_buffer.length() : newpos;
            
            temp = ___global_buffer.substr(pos, newpos - pos);
            
            *(data + i * m + j) = stod(temp);
            
            pos = newpos + 1;
        }
        
        ++i;
    }
    
    return data;
}

vector<double> getNextValue(uint64_t m)
{
    vector<double> newDataInput(m);
    
    if (!___global_file.eof())
    {
        getline(___global_file, ___global_buffer);
        
        size_t pos = 0;
        
        string temp;
        for (uint64_t j = 0; j < m; ++j)
        {
            // TBA
            size_t newpos;
            newpos = ___global_buffer.find(' ', pos + 1);
            
            newpos = newpos == string::npos ? ___global_buffer.length() : newpos;
            
            temp = ___global_buffer.substr(pos, newpos - pos);
            
            newDataInput[j] = stod(temp);
            
            pos = newpos + 1;
        }
    }
    else
    {
        return vector<double>(0);
    }
    
    return newDataInput;
}

void closeFile()
{
    ___global_file.close();
}

//TBA

void writeTime(std::string out, uint64_t n, uint64_t m, int64_t time)
{
    ofstream out_file;
    out_file.open(out, ios::out | ios::app);
    
    out_file << n << "\t" << m << "\t" << time << endl;
    
    out_file.close();
}

void writeMemory(std::string out, uint64_t n, uint64_t m, uint64_t memory)
{
    ofstream out_file;
    out_file.open(out, ios::out | ios::app);
    
    out_file << n << "\t" << m << "\t" << memory << endl;
    
    out_file.close();
}

void writePrecision(std::string out, uint64_t n, uint64_t m, double precision)
{
    ofstream out_file;
    out_file.open(out, ios::out | ios::app);
    
    out_file << "PREC: " << n << "\t" << m << "\t" << precision << endl;
    
    out_file.close();
}

void writeRecovery(std::string &out, uint64_t n, uint64_t m, int64_t result,
                   std::vector<Algorithms::MissingBlock> *missingBlocks)
{
    ofstream out_file;
    out_file.open(out, ios::out | ios::app);
    
    out_file << "Recovery in " << n << " x " << m << " matrix performed in : " << result << "\t" << "ms" << endl
             << "Recovered blocks:" << endl;
    
    for (auto mblock : *missingBlocks)
    {
        Algebra::Vector v = mblock.extractBlock();
        out_file << "MB_" << mblock.column << "," << mblock.startingIndex << "," << mblock.blockSize << " ="
                 << v.toString() << std::endl;
    }
    
    out_file.close();
}

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

MatrixReader::MatrixReader(std::string &input, char sep)
        : file(input),
          separator(sep),
          fileopen(true)
{
    if (!file.is_open())
    {
        std::cout << "Can't open file " << input << std::endl;
        fileopen = false;
    }
}

bool MatrixReader::isValid()
{
    return fileopen;
}

Algebra::Matrix MatrixReader::getFullMatrix()
{
    setFirstRow();
    uint64_t m = rowContainer.size();
    Algebra::Matrix mat(1, m, false);
    
    for (uint64_t j = 0; j < m; ++j)
    {
        mat(0, j) = rowContainer[j];
    }
    
    while (hasNextLine())
    {
        setNextRow();
        mat.append(rowContainer);
    }
    
    return mat;
}

Algebra::Matrix MatrixReader::getFixedMatrix(uint64_t n, uint64_t m)
{
    Algebra::Matrix mat(n, m, false);
    
    rowContainer.reserve(m);
    for (uint64_t j = 0; j < m; ++j)
    {
        rowContainer.push_back(0.0);
    }
    
    for (uint64_t i = 0; i < n; ++i)
    {
        setNextRow();
        for (uint64_t j = 0; j < m; ++j)
        {
            mat(i, j) = rowContainer[j];
        }
    }
    
    return mat;
}

Algebra::Matrix MatrixReader::getFixedRowMatrix(uint64_t n)
{
    setFirstRow();
    uint64_t m = rowContainer.size();
    
    Algebra::Matrix mat(n, m, false);
    
    for (uint64_t j = 0; j < m; ++j)
    {
        mat(0, j) = rowContainer[j];
    }
    
    for (uint64_t i = 1; i < n; ++i)
    {
        setNextRow();
        for (uint64_t j = 0; j < m; ++j)
        {
            mat(i, j) = rowContainer[j];
        }
    }
    
    return mat;
}

Algebra::Matrix MatrixReader::getFixedColumnMatrix(uint64_t m)
{
    Algebra::Matrix mat(1, m, false);
    
    rowContainer.reserve(m);
    for (uint64_t j = 0; j < m; ++j)
    {
        rowContainer.push_back(0.0);
    }
    
    setNextRow(); // in fact first, but we know m, so it's a different call
    for (uint64_t j = 0; j < m; ++j)
    {
        mat(0, j) = rowContainer[j];
    }
    
    while (hasNextLine())
    {
        setNextRow();
        mat.append(rowContainer);
    }
    
    return mat;
}

bool MatrixReader::hasNextLine()
{
    return file.peek() != EOF;
}

Algebra::Vector MatrixReader::readNextLine()
{
    setNextRow();
    
    return Algebra::Vector(rowContainer);
}

void MatrixReader::detectSeparator()
{
    if (buffer.find("\t") != std::string::npos)
    {
        separator = '\t';
    }
    else if (buffer.find(",") != std::string::npos)
    {
        separator = ',';
    }
    else
    {
        separator = ' ';
    }
}

void MatrixReader::setFirstRow()
{
    getline(file, buffer);
    size_t pos = 0;

    if (separator == '\0')
    {
        detectSeparator();
    }
    
    while (pos < buffer.size())
    {
        size_t newpos;
        newpos = buffer.find(separator, pos + 1);
        
        newpos = newpos == string::npos ? buffer.length() : newpos;
        
        std::string temp = buffer.substr(pos, newpos - pos);
        
        if (temp.empty())
        { break; }
        
        rowContainer.push_back(stod(temp));
        
        pos = newpos + 1;
    }
}

void MatrixReader::setNextRow()
{
    uint64_t m = rowContainer.size();
    
    getline(file, buffer);
    size_t pos = 0;

    if (separator == '\0')
    {
        detectSeparator();
    }
    
    for (uint64_t j = 0; j < m; ++j)
    {
        size_t newpos;
        newpos = buffer.find(separator, pos + 1);
        
        newpos = newpos == string::npos ? buffer.length() : newpos;
        
        std::string temp = buffer.substr(pos, newpos - pos);
        
        rowContainer[j] = stod(temp);
        
        pos = newpos + 1;
    }
}

void exportAnyPrecision(std::string &output, uint64_t n, uint64_t m, double precision)
{
    ofstream out_file;
    out_file.open(output, ios::out | ios::app);
    
    out_file << n << "\t" << m << "\t" << precision << endl;
    
    out_file.close();
}

void exportAnyRuntime(std::string &output, uint64_t n, uint64_t m, int64_t runtime)
{
    ofstream out_file;
    out_file.open(output, ios::out | ios::app);
    
    out_file << n << "\t" << m << "\t" << runtime << endl;
    
    out_file.close();
}

void exportDecompOutput(std::string &output, const Algebra::Matrix &Load, const Algebra::Matrix &Rel,
                        const std::vector<double> &centroidValues)
{
    exportMatrix(output + "_L.txt", Load);
    
    exportMatrix(output + "_R.txt", Rel);
    
    ofstream out_file;
    // Centroid Values
    out_file.open(output + "_centroid.txt", ios::out);
    for (double elem : centroidValues)
    {
        out_file << elem << std::endl;
    }
    out_file.close();
}

void exportMatrix(std::string output, const Algebra::Matrix &mx)
{
    ofstream out_file;
    out_file.open(output, ios::out);
    
    for (uint64_t i = 0; i < mx.dimN(); ++i)
    {
        for (uint64_t j = 0; j < mx.dimM() - 1; ++j)
        {
            out_file << mx(i, j) << " ";
        }
        out_file << mx(i, mx.dimM() - 1) << std::endl;
    }
    out_file.close();
}

} // namespace MathIO
