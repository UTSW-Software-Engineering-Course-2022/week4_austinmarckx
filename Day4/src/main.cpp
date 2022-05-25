
/**
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <execution>

 // TODO: support multiple chromosome
class FReference
{
public:
    // Default constructor
    FReference() {};

    // Constructor with fasta filename
    FReference(const std::string& InFilename) {
        LoadFromFasta(InFilename);
    }
    // Constructor with fasta filename and suffixarray file
    FReference(const std::string& InFilename, const std::string& saFilename) {
        LoadFromFasta(InFilename);
        LoadSuffixArray(saFilename);
    }

    void LoadSuffixArray(const std::string& saFilename) {
        std::ifstream inputFile(saFilename);
        std::vector<int> suffArr;
        int element;
        while (inputFile >> element)
        {
            suffArr.push_back(element);
        }
        SetSuffixArray(suffArr, 0);
    }

    void LoadFromString(const std::string& InSequence, int seqIndex = 0)
    {
        Sequence[seqIndex] = InSequence;
        if (Sequence[seqIndex].back() != '$') {
            Sequence[seqIndex].push_back('$');
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
        SequenceName.clear();

        int seqNum = -1;
        while (getline(fs, buf)) {
            if (buf.length() == 0) {
                // skip empty line
                continue;
            }
            // If Find a sequence, make space and get name
            if (buf[0] == '>') {
                // header line
                seqNum++;
                SequenceName.push_back(buf.substr(1, buf.length() - 1));
                Sequence.push_back("");
                SequenceLength.push_back(0);
                suffixArray.push_back({ 0 });
            }
            // Read in sequence
            else
                Sequence[seqNum].append(buf);
        }
        // append '$' as End of Sequence mark
        for (int i = 0; i <= seqNum; i++) {
            Sequence[i].append("$");
            SequenceLength[i] = Sequence[i].length();
        }
    }

    // Save to file
    void SaveIndexesToFile(const std::string& OutFilename)
    {
        // initialize output file
        std::ofstream outFile;
        outFile.open(OutFilename);

        for (int seqNum = 0; seqNum < Sequence.size(); seqNum++)
            for (int i = 0; i < Sequence[seqNum].size(); i++) {
                // Current suffex has index i+1, and is the substring from i to end.
                outFile << suffixArray[seqNum][i] << std::endl;
            }
        // close file
        outFile.close();
    }

    // Get Seq Len
    size_t GetSequenceLength(int seqIndex) const {
        return SequenceLength[seqIndex];
    }
    // Set suffix array
    void SetSuffixArray(std::vector<int> inputSuffArray, int seqIndex) {

        suffixArray[seqIndex] = inputSuffArray;
        std::cout << "Suffix Array Set" << std::endl;
    }

public:
    std::vector<std::string> SequenceName;
    std::vector<std::string> Sequence;
    std::vector<size_t> SequenceLength;
    std::vector<std::vector<int>> suffixArray;
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

// Struct to store suffixes
struct suffix
{
    // inspiration: https://github.com/saadtaame/suffix-sort
    // Also: https://github.com/UTSW-Software-Engineering-Course-2022/week_4_aleks
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
bool suffixComparision(const suffix& Apples, const suffix& Oranges) {

    if (Apples.rank[0] < Oranges.rank[0])
        return true;

    else if (Apples.rank[0] == Oranges.rank[0])
        return (Apples.rank[1] < Oranges.rank[1]);

    return false;
}

/* Create Suffix Array from Sequence

*/
std::vector<int> SuffixArrayFromSequence(FReference refSeq, int refSeqIndex = 0)
{
    // Get length
    int strLen = (int)refSeq.Sequence[refSeqIndex].size();

    // Instantiate output suffix array; fill with original indexes
    std::vector<int> rankArray(strLen);
    for (int i = 0; i < strLen; i++) {
        rankArray[i] = (int)refSeq.Sequence[refSeqIndex][i];
    }

    // Instantiate Suffix vector, 
    std::vector<suffix> suffVect(strLen);

    // Go over suffix vector and sort ranks
    for (int len = 1; ; len *= 2) {
        for (int j = 0; j < strLen; j++) {
            suffVect[j] = suffix(
                j,
                rankArray[j],
                ((j + len < strLen) ? rankArray[j + len] : -1)
            );
        }
        // Sort
        sort(suffVect.begin(), suffVect.end(), suffixComparision);

        // Update suffix array from sorted suffix vector
        rankArray[suffVect[0].index] = 0;
        for (int i = 1; i < strLen; i++) {
            rankArray[suffVect[i].index] = rankArray[suffVect[i - 1].index] + (suffixComparision(suffVect[i - 1], suffVect[i]) ? 1 : 0);
        }

        // Early exit condition
        if (rankArray[suffVect.back().index] == strLen - 1)
            break;
    }
    // Extract suffix array output from sorted suffix vector
    std::vector<int> suffixArray(strLen);
    for (int i = 0; i < strLen; i++) {
        suffixArray[i] = suffVect[i].index;
    }
    std::cout << "Suffix Array determined..." << std::endl;
    // Return
    return suffixArray;
}

// struct to hold alignments
struct alignment {
    std::string querySeqName;
    std::string refSeqName;
    int refSeqIndex;
    // Default constructor
    alignment() {}

    // constructor
    alignment(std::string qName, std::string rName, int rIndex) {
        querySeqName = qName;
        refSeqName = rName;
        refSeqIndex = rIndex;
    }
};


void SaveAlignmentsToFile(const std::string& OutFilename, const std::vector<alignment> aln)
{
    // initialize output file
    std::ofstream outFile;
    outFile.open(OutFilename);
    for (int i = 0; i < aln.size(); i++) {
        // Current suffex has index i+1, and is the substring from i to end.
        outFile << aln[i].querySeqName << "\t" <<
            aln[i].refSeqName << "\t" << aln[i].refSeqIndex << std::endl;
    }
    // close file
    outFile.close();
}

void PrintAlignmentsToTerminal(const std::vector<alignment> aln)
{
    // initialize output file
    for (int i = 0; i < aln.size(); i++) {
        // Current suffex has index i+1, and is the substring from i to end.
        std::cout << aln[i].querySeqName << "\t" <<
            aln[i].refSeqName << "\t" << aln[i].refSeqIndex << std::endl;
    }
}

// Midpoint for Binary search
int GetMidpoint(int left, int right) {
    return(int((left + right) / 2));
}

std::vector<alignment> AlignQueryToSuffixArray(FReference ref, FReference query, int querySeqIndex, int refSeqIndex = 0) {
    // Setup
    int queryLength = query.GetSequenceLength(querySeqIndex) - 1;
    int refSeqLength = query.GetSequenceLength(refSeqIndex) - 1;
    int left = 0;
    int right = ref.GetSequenceLength(refSeqIndex);
    int mid = GetMidpoint(left, right);
    std::vector<int> suffArrayCpy = ref.suffixArray[refSeqIndex];
    std::string querySeqCpy = query.Sequence[querySeqIndex].substr(0, queryLength);

    std::cout << "querySeqCpy: " << querySeqCpy << std::endl;

    // output
    std::vector<alignment> matches = {};
    std::string tmp;

    // Try and use seed of query to increase comparison speed
    int seedSize = 1;
    std::string querySeqSubStr = querySeqCpy.substr(0, seedSize);
    
    // Binary search
    while (left != right) {
        // get string at midpoint of suffix array
        tmp = ref.Sequence[refSeqIndex].substr(suffArrayCpy[mid], seedSize);
        std::cout << "tmp: " << tmp << "  mid: " << mid << "  r: " << right << "  l: " << left << std::endl;

        // Too low, move left pointer
        if (tmp < querySeqSubStr) {
            if (left == mid) {
                left++;
            }
            else {
                left = mid;
            }

        }// Too high, move right pointer 
        else if (tmp > querySeqSubStr) {
            if (right == mid) {
                right--;
            }
            else {
                right = mid;
            }
        } // O/W match
        else {
            // increase seed Size
            if (seedSize != queryLength) {
                seedSize *= 2;
                // If seed size too big, set to query length
                if (seedSize > queryLength) {
                    seedSize = queryLength;
                }
                // update query substring
                querySeqSubStr = querySeqCpy.substr(0, seedSize);
            }
            else { // if match on max seed size: add match
                // Add the entry to output
                matches.push_back(alignment(query.SequenceName[querySeqIndex], ref.SequenceName[refSeqIndex], suffArrayCpy[mid] + 1)); // 0->1 Index
                // remove the match from suffixArray
                suffArrayCpy.erase(suffArrayCpy.begin() + mid);
            }
        }

        // recalculate midpoint
        mid = GetMidpoint(left, right);
    }

    return matches;
}


struct BWTArray {
    std::string EOS;
    int A;
    int C;
    int G;
    int T;
    std::vector<int> rank;
    std::string first;
    std::string last;

    BWTArray() {
        std::string EOS = "$";
        int A = 0;
        int C = 0;
        int G = 0;
        int T = 0;
        std::vector<int> rank; //character, counts 
        std::string first = "";
        std::string last = "";
    }
};

BWTArray BWTFromSuffixArray(FReference seq, int seqIndex) {
    BWTArray output = BWTArray();

    int seqLen = seq.suffixArray[seqIndex].size();

    for (int i = 0; i < seqLen; i++) {
        if (seq.suffixArray[seqIndex][i] == 0) {
            output.last.append("$");
        }
        else {
            output.last.append(seq.Sequence[seqIndex].substr(seq.suffixArray[seqIndex][i] - 1, 1));
        }
    }
    return output;
}

void CalculateBWTRank(BWTArray &bwt) {
    // Calculate rank
    bwt.rank.resize(bwt.last.length());
    char matches[5] = { '$','A','C','G','T'};
    std::vector<int> counts(5);
    std::vector<int> totalCounts(5);
    std::vector<int> rankIndex(5);

    // Count and store
    bwt.A = std::count(bwt.last.begin(), bwt.last.end(), 'A');
    bwt.C = std::count(bwt.last.begin(), bwt.last.end(), 'C');
    bwt.G = std::count(bwt.last.begin(), bwt.last.end(), 'G');
    bwt.T = std::count(bwt.last.begin(), bwt.last.end(), 'T');
    totalCounts[0] = 1;
    totalCounts[1] = bwt.A;
    totalCounts[2] = bwt.C;
    totalCounts[3] = bwt.G;
    totalCounts[4] = bwt.T;

    rankIndex[0] = bwt.last.length()-totalCounts[0];
    rankIndex[4] = rankIndex[0] - totalCounts[4];
    rankIndex[3] = rankIndex[4] - totalCounts[3];
    rankIndex[2] = rankIndex[3] - totalCounts[2];
    rankIndex[1] = rankIndex[2] - totalCounts[1];

    for (int i = 0; i < bwt.last.length(); i++) {
        // Get current char
        char tmpChar = bwt.last[i];
        // find matching char, update counts and rank
        for(int j = 0; j < 5; j++)
            if (tmpChar == matches[j]) {
                counts[j]++;
                bwt.rank[i] = rankIndex[j] + counts[j];
                break;
            }
        //std::cout << "$: " << counts[0] << " A: " << counts[1] << " C: " << counts[2] << " G: " << counts[3] << " T: " << counts[4] << std::endl;
    }
}

std::string SequenceFromBWTAndSuffixArray(FReference seq, int seqIndex, BWTArray bwt) {
    // Calculate rank

    std::cout << "Calculating Rank for BWT..." << std::endl;
    CalculateBWTRank(bwt);
    std::cout << "Rank for BWT Calculated." << std::endl;
    
    // create output string 
    std::string output = "";
    int r = 0;
    char c = bwt.last[r];

    while (c != '$') {
        //std::cout << "r: " << r << " c: " << c << std::endl;
        output += c;
        //std::cout << "rank " << bwt.rank[r] << " next c: " << bwt.last[r] << std::endl;
        //std::cout << "rank " << bwt.rank[58] << " next c: " << bwt.last[58] << std::endl;

        r = bwt.rank[r];
        c = bwt.last[r];
    }
    //std::cout << "Output generated... " << output << " " << output.length() << std::endl;
    
    // For some reason this generates the reverse of the input sequence... unreverse it.
    std::reverse(output.begin(), output.end());
    return output += '$';
}


int main(int argc, char* argv[])
{

    // sa InReferenceFastaFile OutSuffixArrayFile
    if (argc < 3) {
        PrintUsage(GetFilename(argv[0]));
        return 1;
    }

    {
        // Load reference and suffix array
        std::cout << "Loading Reference Suffix Array..." << std::endl;
        FReference ref;
        ref.suffixArray.push_back({ 0 });
        ref.LoadSuffixArray(argv[1]);
        std::cout << "Reference Suffix Array Loaded." << std::endl;

        // Load BWT
        std::cout << "Loading BWT..." << std::endl;
        BWTArray bwt;
        std::ifstream read;
        read.open(argv[2]);
        std::getline(read, bwt.last);
        read.close();
        std::cout << "BWT Loaded." << std::endl;

        // Recover original sequence from BWT and suffix Array and time it
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::cout << "Recovering Sequence from BWT and Suffix Array..." << std::endl;

        ref.Sequence.push_back({ 0 });
        ref.Sequence[0] = SequenceFromBWTAndSuffixArray(ref, 0, bwt);
        ref.SequenceLength.push_back(0);
        ref.SequenceLength[0] = ref.Sequence[0].length();
        ref.SequenceName.push_back("");
        ref.SequenceName[0] = "smalls";

        // Print and time output 
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto durationFPextraction = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Sequence recovery time: " << durationFPextraction << std::endl;

        // Load Queries
        std::cout << "Loading Queries..." << std::endl;
        FReference query(argv[3]);
        std::cout << "Queries Loaded." << std::endl;

        std::cout << query.Sequence[0] << std::endl;
        std::cout << ref.Sequence[0] << std::endl;


        std::cout << "Aligning Queries..." << std::endl;
        std::vector<alignment> alignments;

        // All queries
        for (int i = 0; i < query.Sequence.size(); i++) {
            // 1st 1000 queries
            //for (int i = 0; i <1000; i++) {
                // Calc Time step per loop: start
            std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
            //std::cout << "attempt to align";
            std::vector<alignment> tmpAlignments = AlignQueryToSuffixArray(ref, query, i);
            // std::cout << "got temp align";
            // add new matches to alignment vector
            for (int j = 0; j < tmpAlignments.size(); j++) {
                alignments.push_back(tmpAlignments[j]);
            }

            // Time step per iteration: end 
            std::cout << "Alignment #" << i + 1 << std::endl;
            std::cout << "Num Alignments Found: " << tmpAlignments.size() << std::endl;

            std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
            auto durationFPextraction = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
            std::cout << "Alignment time: " << durationFPextraction << std::endl;

        }
        std::cout << "Queries Aligned..." << std::endl;

        // Print alignments to file
        //PrintAlignmentsToTerminal(alignments);

        // Save indexes to output file
        std::cout << "Saving alignments..." << std::endl;
        SaveAlignmentsToFile(GetFilename(argv[4]), alignments);
        std::cout << "Alignments saved." << std::endl;
        
        std::cout << "Program completed." << std::endl;
    }

    return 0;
}

