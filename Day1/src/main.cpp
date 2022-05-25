/**
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <regex>

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
        // initialize output file
        std::ofstream outFile;
        outFile.open(OutFilename);
        for (int i = 0; i < Sequence.length(); i++) {
            // Current suffex has index i+1, and is the substring from i to end.
            outFile << suffixArray[i] << std::endl;
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
    void SetSuffixArray(std::vector<int> inputSuffArray) {
        suffixArray = inputSuffArray;
    }
 
public:
    std::string Sequence;
    size_t SequenceLength;
    std::vector<int> suffixArray;
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

// inspiration: https://github.com/saadtaame/suffix-sort
// Also: https://github.com/UTSW-Software-Engineering-Course-2022/week_4_aleks
// Struct to store suffixes
struct suffix
{
    int index;
    int rank[2];

    // suffix constructor
    suffix() {}
    suffix(int index_, int rank_one, int rank_two) {
        index = index_;
        rank[0] = rank_one;
        rank[1] = rank_two;
    }
};

// comparision for sort:
bool suffixComparision(const suffix &Apples, const suffix &Oranges) {
    
    if (Apples.rank[0] < Oranges.rank[0])
        return true;

    else if (Apples.rank[0] == Oranges.rank[0])
        return (Apples.rank[1] < Oranges.rank[1]);

    return false;
}

/* Create Suffix Array from Sequence

*/
std::vector<int> SuffixArrayFromSequence(FReference refSeq)
{
    // Get length
    int strLen = (int)refSeq.Sequence.size();
    
    //std::cout << strLen << std::endl;

    // Instantiate output suffix array; fill with original indexes
    std::vector<int> rankArray(strLen);
    for (int i = 0; i < strLen; i++) {
        rankArray[i] = (int)refSeq.Sequence[i];
        
        //std::cout << rankArray[i] <<std::endl;
    }
    //std::cout << "Passed Rank Array Instantiation" << std::endl;

    // Instantiate Suffix vector, 
    std::vector<suffix> suffVect(strLen);

    // Go over suffix vector and sort ranks
    for (int len = 1; ; len *= 2) {
        for (int j = 0; j < strLen; j++) {
            suffVect[j] = suffix(
                j,
                rankArray[j],
                ((j + len < strLen) ? rankArray[j+len] : -1)
            );

            //std::cout << suffVect[j].index << " " << suffVect[j].rank[0] << " " << suffVect[j].rank[1] << std::endl;
        }
        //std::cout << "Passed SuffixVector Instantiation" << std::endl;

        sort(suffVect.begin(), suffVect.end(), suffixComparision);

        //std::cout << "Passed Sort" << std::endl;
        // Update suffix array from sorted suffix vector
        
        // THERE IS A BUG HERE
        //std::cout << suffVect[0].index << std::endl;

        rankArray[suffVect[0].index] = 0;
        //std::cout << "Enter rank update loop" << std::endl;
        for (int i = 1; i < strLen; i++) {
            rankArray[suffVect[i].index] = rankArray[suffVect[i - 1].index] + (suffixComparision(suffVect[i - 1], suffVect[i]) ? 1 : 0);
        }
        //std::cout << "Passed Suffix Array Update" << std::endl;

        // Early exit condition
        if (rankArray[suffVect.back().index] == strLen - 1)
            break;
    } 
    //std::cout << "Passed ALL sort loop" << std::endl;

    std::vector<int> suffixArray(strLen);
    for (int i = 0; i < strLen; i++) {
        suffixArray[i] = suffVect[i].index;
    }
    //std::cout << "Suffix Array Created" << std::endl;
    // Return
    return suffixArray;
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
        ref.SetSuffixArray(SuffixArrayFromSequence(ref));
        std::cout << "Suffix array created..." << std::endl;

        // Save indexes to output file
        std::cout << "Saving indexes..." << std::endl;
        ref.SaveIndexesToFile(GetFilename(argv[2]));
        std::cout << "Indexes Saved to "<< GetFilename(argv[2]) << std::endl;
        std::cout << "Program completed." << std::endl;
    }

    return 0;
}
