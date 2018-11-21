//
// Created by Zakhar on 03.04.2017.
//

#pragma once

#include <string>
#include <iostream>

using namespace std;

void printUsage()
{
    cout << endl << "Usage: ./incCD -obligatory arg [-optional arg]" << endl
         << endl
         << "-action {arg}, -act {arg}" << endl
         << "    | arg:" << endl
         << "        | dec, d  - centroid decomposition" << endl
         << "        | rec, r  - recovery of missing values" << endl
         << "        | norm, n - normalization of the matrix" << endl
         << "    | the baseline action to perform" << endl
         << endl
         << "-test {arg}, -t {arg}" << endl
         << "    | arg:" << endl
         << "        | out, o      - the result of the <action>" << endl
         << "        | prec, p     - precision of the <action>" << endl
         << "        | runtime, rt - runtime of the <action>" << endl
         << "    | choose what to output from the action" << endl
         << endl
         << "-input {str}, -in {str}" << endl
         << "    | file name, where to take the input matrix from" << endl
         << endl
         << "-output {str}, -out {str}" << endl
         << "    | file name, where to store the result of <test>" << endl
         << endl
         << "[-n {int}] default(0)" << endl
         << "    | amount of rows to load from the input file" << endl
         << "    | 0 - load all of them" << endl
         << endl
         << "[-m {int}] default(0)" << endl
         << "    | amount of columns to load from the input file" << endl
         << "    | 0 - load all of them" << endl
         << endl
         << "[-k {int}] default(m)" << endl
         << "    | amount of columns of truncated decomposition to keep" << endl
         << "    | 0 (dec) - will be set to be equal to m" << endl
         << "    | 0 (rec) - will be automatically detected" << endl
         << endl
         << "[-sv {arg}] default(lsv-base)" << endl
         << "    | arg:" << endl
         << "        | issv       - incremental scalable sign vector" << endl
         << "        | issv+      - ISSV+ (lifted condition of 1 update per iteration)" << endl
         << "        | issv-init  - ISSV  but with linear init" << endl
         << "        | issv+-init - ISSV+ but with linear init" << endl
         << "        | lsv-base   - local sign vector" << endl
         << "        | lsv-noinit - LSV with skipped init" << endl
         << "    | a sign vector search algorithm to be used with CD, applicable for <action> dec & rec" << endl
         << "    | invalid value will result in a warning message and a fallback to a default algorithm" << endl
         << endl
         << "[-xtra {arg}] default(<none>)" << endl
         << "    | arg: " << endl
         << "        | inc, i         - turn on incremental mode (supplying -n and -imax is obligatory)" << endl
         << "        | ilastonly, ilo - (incremental) print only the last action of the incremental <action>" << endl
         << "        | full           - (non-incremental) export rt/prec with info about n, m" << endl
         << "        | batch, bat     - (recovery) disable sign vector caching during recovery" << endl
         << "        | recnorm        - (recovery) use z-score normalization during recovery" << endl
         << "    | turns on a special flag" << endl
         << "    | multiple flags require \"-xtra flag1 -xtra flag2 ...\"" << endl
         << endl
         << "[-istep {int}] default(1)" << endl
         << "    | (inc) amount of rows to load before re-applying <action>" << endl
         << endl
         << "[-imax {int}] default(0)" << endl
         << "    | (inc) maximum amount of rows to apply <action> to" << endl
         << "    |       0 - load as much as possible" << endl
         << endl
         << "[-rcdopt {int}] default(0)" << endl
         << "    | (rec) apply optimization with the resp. code" << endl
         << "    |       overrides \"-xtra bat\"" << endl
         << endl;
}

enum class PAction
{
    Decomposition, Recovery, Normalization, Undefined
};

enum class PTestType
{
    Output, Precision, Runtime, Undefined
};

struct PTestParam
{
    bool Incremental, IncrementalLastOnly, UseBatchCD, UseFullOutput, RecoverNormalized;
    
    explicit PTestParam()
            : Incremental(false), IncrementalLastOnly(false),
              UseBatchCD(false), UseFullOutput(false), RecoverNormalized(false)
    { }
};

int CommandLine2(
        int argc, char *argv[],
        
        PAction &action, PTestType &test,
        PTestParam &testParams,
        
        std::string &input, std::string &output,
        
        uint64_t &n, uint64_t &m, uint64_t &k,
        uint64_t &inc, uint64_t &max,
        
        uint64_t &rcdopt, string &signVector
)
{
    string temp;
    for (int i = 1; i < argc; ++i)
    {
        temp = argv[i];
        
        // man request
        if (temp.empty())
        {
            continue;
        }
        else if (temp == "--help" || temp == "-help" || temp == "/?")
        {
            printUsage();
            return 1;
        }
            
            // Action type, base algorithm to be used
        else if (temp == "-act" || temp == "-action")
        {
            ++i;
            temp = argv[i];
            
            if (temp == "dec" || temp == "d")
            {
                action = PAction::Decomposition;
            }
            else if (temp == "rec" || temp == "r")
            {
                action = PAction::Recovery;
            }
            else if (temp == "norm" || temp == "n")
            {
                action = PAction::Normalization;
            }
            else
            {
                cout << "Unrecognized -action argument" << endl;
                printUsage();
                return EXIT_FAILURE;
            }
        }
            
            // Test type, how the algorithm is applied
        else if (temp == "-test" || temp == "-t")
        {
            ++i;
            temp = argv[i];
            
            if (temp == "out" || temp == "o")
            {
                test = PTestType::Output;
            }
            else if (temp == "prec" || temp == "p")
            {
                test = PTestType::Precision;
            }
            else if (temp == "runtime" || temp == "rt")
            {
                test = PTestType::Runtime;
            }
            else
            {
                cout << "Unrecognized -test argument" << endl;
                printUsage();
                return EXIT_FAILURE;
            }
        }
            
            // Additional parameters
        else if (temp == "-sv")
        {
            ++i;
            signVector = argv[i];
        }
        
        else if (temp == "-xtra")
        {
            ++i;
            temp = argv[i];
            
            if (temp == "inc" || temp == "i")
            {
                testParams.Incremental = true;
            }
            else if (temp == "batch" || temp == "bat")
            {
                testParams.UseBatchCD = true;
            }
            else if (temp == "full")
            {
                testParams.UseFullOutput = true;
            }
            else if (temp == "ilastonly" || temp == "ilo")
            {
                testParams.IncrementalLastOnly = true;
            }
            else if (temp == "recnorm")
            {
                testParams.RecoverNormalized = true;
            }
            else
            {
                cout << "Unrecognized -xtra argument" << endl;
                printUsage();
                return EXIT_FAILURE;
            }
        }
            
            // in/out
        else if (temp == "-input" || temp == "-in")
        {
            ++i;
            input = argv[i];
        }
        else if (temp == "-output" || temp == "-out")
        {
            ++i;
            output = argv[i];
        }
            
            // Dimensions to override defaults
        else if (temp == "-n")
        {
            ++i;
            temp = argv[i];
            
            n = static_cast<uint64_t>(stoll(temp));
        }
        else if (temp == "-m")
        {
            ++i;
            temp = argv[i];
            
            m = static_cast<uint64_t>(stoll(temp));
        }
        else if (temp == "-k")
        {
            ++i;
            temp = argv[i];
            
            k = static_cast<uint64_t>(stoll(temp));
        }
            
            // Incremental options
        else if (temp == "-istep")
        {
            ++i;
            temp = argv[i];
            
            inc = static_cast<uint64_t>(stoll(temp));
        }
        else if (temp == "-imax")
        {
            ++i;
            temp = argv[i];
            
            max = static_cast<uint64_t>(stoll(temp));
        }
            
            // Other values
        else if (temp == "-rcdopt")
        {
            ++i;
            temp = argv[i];
            
            rcdopt = static_cast<uint64_t>(stoll(temp));
        }
        
        else
        {
            cout << "Unrecognized CLI parameter" << endl;
            printUsage();
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}