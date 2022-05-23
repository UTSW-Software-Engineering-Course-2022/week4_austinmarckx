/**
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

// Struct to store suffixes
struct suffix
{
    int index;
    std::string suff;
};

// TODO: support multiple chromosome
class FReference 
{
public:
    // Default constructor
    FReference() {};

    // Constructor with fasta filename
    FReference(const std::string& InFilename) {
        LoadFromFasta(InFilename);
        SetSequenceLength();
    }

    void LoadFromString(const std::string& InSequence)
    {
        Sequence = InSequence;
        if (Sequence.back() != '$') {
            Sequence.push_back('$');
        }
    }
    
    void LoadFromFasta(const std::string& InFilename)
    {
        std::ifstream fs(InFilename);
        if (!fs) {
            std::cerr << "Can't open file: " << InFilename << std::endl;
            return;
        }

        std::string buf;

        Sequence.clear();

        while(getline(fs, buf)) {
            if (buf.length() == 0) {
                // skip empty line
                continue;
            }

            if (buf[0] == '>') {
                // header line
                // TODO: save chromosome name
                continue;
            }

            Sequence.append(buf);
        }

        // append '$' as End of Sequence mark
        Sequence.append("$");

    }

    // Save to file
    void SaveIndexesToFile(const std::string& OutFilename)
    {
        // Get indexes
        std::vector<int> indexes = GetSuffixArrayIndex();
        // initialize output file
        std::ofstream outFile;
        outFile.open(OutFilename);
        for (int i = 0; i < indexes.size(); i++) {
            // Current suffex has index i+1, and is the substring from i to end.
            outFile << indexes[i] << std::endl;
        }
        // close file
        outFile.close();
    }

    // Get and Set sequence length
    void SetSequenceLength() {
        SequenceLength = Sequence.length();
    }
    size_t GetSequenceLength() const {
        return SequenceLength;
    } 
    // Set suffix array
    void SetSuffixArray() {
        suffixArray = SuffixArrayFromSequence();
    }
    // return vector of ints continain the suffix array indexes
    std::vector<int> GetSuffixArrayIndex() {
        std::vector<int> indexes;
        indexes.resize(suffixArray.size());
        for (int i = 0; i < suffixArray.size(); i++) {
            indexes[i] = suffixArray[i].index;
        }
        return indexes;
    }

    // Comparision for sort:
    static int suffixComparision(suffix one, suffix two) {
        return std::strcmp(one.suff.c_str(), two.suff.c_str()) < 0 ? 1 : 0;
    }

    /* Create Suffix Array from Sequence
    
    */
    std::vector<suffix> SuffixArrayFromSequence()
    {   
        // Allocate for suffix array
        suffixArray.resize(SequenceLength);

        // Fill Suffix array
        for (int i = 0; i < SequenceLength; i++) {
            // Current suffex has index i+1, and is the substring from i to end.
            suffixArray[i].index = i;
            suffixArray[i].suff = Sequence.substr(i, SequenceLength);
        }
        // Sort
        std::sort(suffixArray.begin(), suffixArray.end(), suffixComparision);

        // Return
        return suffixArray;
    }

public:
    std::string Sequence;
    size_t SequenceLength;
    std::vector<suffix> suffixArray;
};


/* 
 * Return filename from full path of file.
 *
 * InFilename   [In] Full path of file
 *
 * Return
 *   base filename
 *
 */
static std::string GetFilename(const std::string& InFilename)
{
    const size_t pos = InFilename.find_last_of("/\\");
    if (pos == std::string::npos) {
        return InFilename;
    }

    return InFilename.substr(pos + 1);
}

void PrintUsage(const std::string& InProgramName)
{
    std::cerr << "Invalid Parameters" << std::endl;
    std::cerr << "  " << InProgramName << " Reference_Fasta_File SuffixArray_File" << std::endl;
}

/**
 * 
 *
 */
int main(int argc, char* argv[])
{

    // sa InReferenceFastaFile OutSuffixArrayFile
    if (argc < 3) {
        PrintUsage(GetFilename(argv[0]));
        return 1;
    }

    {
        FReference ref(argv[1]);

        std::cout << "Reference sequence length: " << ref.Sequence.length() << std::endl;

        // print first 100bp
        std::cout << ref.Sequence.substr(0, 100) << std::endl;
        ref.SuffixArrayFromSequence();
        std::cout << "Suffix array created..." << std::endl;

        // Save indexes to output file
        std::cout << "Saving indexes..." << std::endl;
        ref.SaveIndexesToFile(GetFilename(argv[2]));
        std::cout << "Indexes Saved to "<< GetFilename(argv[2]) << std::endl;
        std::cout << "Program completed." << std::endl;
    }

    return 0;
}
