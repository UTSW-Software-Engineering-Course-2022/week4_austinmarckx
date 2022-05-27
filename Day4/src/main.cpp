
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
#include <execution>]
#include <set>
#include <iterator>

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
    void SetSuffixArray(std::vector<int> &inputSuffArray, int seqIndex) {

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
std::vector<int> SuffixArrayFromSequence(FReference &refSeq, int refSeqIndex = 0)
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

// Recursive binary suffix array search
std::vector<int> AlignQueryToSeq_SuffixArray(FReference &ref, FReference &query, int left, int right, int querySeqLength, int querySeqIndex, int refSeqIndex = 0) {
    // Setup
    int mid = left + (right - left) / 2;

    // output
    std::vector<int> matches = {};
    matches.reserve(3);
    int tmp;

    // Binary search
    while (left <= right) {       

        // compare query to reconstructed sequence here... ------->                                                                          <--------------------
        tmp = strncmp(query.Sequence[querySeqIndex].c_str(), ref.Sequence[refSeqIndex].c_str() + ref.suffixArray[refSeqIndex][mid], querySeqLength);
        //std::cout << "tmp: " << tmp << "  l: " << left << "  mid: " << mid << "  r: " << right<< std::endl;

        // Too low, move left pointer
        if (tmp > 0) {
            if (left == mid) {
                left++;
            }
            else {
                left = mid + 1;
            }

        }// Too high, move right pointer 
        else if (tmp < 0) {
            if (right == mid) {
                right--;
            }
            else {
                right = mid - 1;
            }
        } // O/W match
        else {
            // Add the entry to output
            //std::cout << "Match Found at: " << mid << std::endl;
            matches.push_back(ref.suffixArray[refSeqIndex][mid] + 1); // 0->1 Index
            
            if (left == mid) {
                left++;
            }
            if (right == mid) {
                right--;
            }
            // Recursively look on left and right of mid:
            //std::cout << "Enter Left Recurse" << std::endl;
            std::vector<int> leftRecurs = AlignQueryToSeq_SuffixArray(ref, query, left, mid - 1, querySeqLength, querySeqIndex);

            //std::cout << "Enter Right Recurse" << std::endl;
            std::vector<int> rightRecurs = AlignQueryToSeq_SuffixArray(ref, query, mid + 1, right, querySeqLength, querySeqIndex);

            // Add left and right recursive to matches
            for (int i = 0; i < leftRecurs.size(); i++) {
                matches.push_back(leftRecurs[i]);
            }
            for (int i = 0; i < rightRecurs.size(); i++) {
                matches.push_back(rightRecurs[i]);
            }
            
            return matches;
        }
        
        // recalculate midpoint
        mid = left + (right - left) / 2;
    }

    return matches;
}


struct BWTArray {
    std::string EOS;
    std::vector<int> rank;
    std::set<char> uniqueChars;
    std::vector<int> rankIndex;
    std::vector<int> counts;
    std::vector<int> totalCounts;
    std::string first;
    std::string last;

    BWTArray() {
        std::string EOS = "$";
        std::set<char> uniqueChars;
        std::vector<int> rank; //character, counts 
        std::vector<int> rankIndex;
        std::string first = "";
        std::string last = "";
        std::vector<int> counts;
        std::vector<int> totalCounts;
    }
};


std::string PartialSequenceFromBWTAndSuffixArray(FReference& seq, int saIndex, BWTArray& bwt, int lenOfRecon) {
    // create output string 
    std::string output = "";
    output.reserve(bwt.last.length());

    // All I have to do is fix this and I think I got it...

    // Find first 

    // saIndex == mid  == the "middle" index of the Suffix Array
    int r = bwt.rank[saIndex];//.last[saIndex];//find(seq.suffixArray[0].begin(), seq.suffixArray[0].end(), saIndex) - seq.suffixArray[0].begin() ;

    // init
    output.append(bwt.last.substr(r, 1));

    for (int i = 0; i < lenOfRecon; i++) {
        //std::cout << "r: " << r << " c: " << bwt.last.substr(bwt.rank[r], 1).c_str() << std::endl; // Debug
        output.append(bwt.last.substr(bwt.rank[r], 1));
        r = bwt.rank[r];
        if (strncmp(bwt.last.substr(bwt.rank[r], 1).c_str(), "$", 1) == 0) {
            break;
        }
    }
    // Fix missing character
    std::reverse(output.begin(), output.end());
    output.erase(0, 1);
    output.append("$");
    return output;
}


// Recursive binary suffix array search
std::vector<int> AlignQueryToBWT_SuffixArray(FReference& ref, FReference& query, BWTArray &bwt, int left, int right, int querySeqLength, int querySeqIndex, int refSeqIndex = 0) {
    // Setup
    int mid = left + (right - left) / 2;

    // output
    std::vector<int> matches = {};
    matches.reserve(3);
    int tmp;

    // Binary search
    while (left <= right) {
        // compare query to reconstructed sequence here... ------->                                                     <--------------------
        tmp = strncmp(query.Sequence[querySeqIndex].c_str(), PartialSequenceFromBWTAndSuffixArray(ref, mid, bwt, querySeqLength).c_str(), querySeqLength);
        std::cout << "tmp: " << tmp << "  l: " << left << "  mid: " << mid << "  r: " << right<< std::endl;
        std::cout << "q: " << query.Sequence[querySeqIndex].c_str() << "  ref: " << PartialSequenceFromBWTAndSuffixArray(ref, mid, bwt, querySeqLength).c_str() << std::endl;

        // Too low, move left pointer
        if (tmp > 0) {
            if (left == mid) {
                left++;
            }
            else {
                left = mid + 1;
            }

        }// Too high, move right pointer 
        else if (tmp < 0) {
            if (right == mid) {
                right--;
            }
            else {
                right = mid - 1;
            }
        } // O/W match
        else {
            // Add the entry to output
            //std::cout << "Match Found at: " << mid << std::endl;
            matches.push_back(ref.suffixArray[refSeqIndex][mid] + 1); // 0->1 Index

            if (left == mid) {
                left++;
            }
            if (right == mid) {
                right--;
            }
            // Recursively look on left and right of mid:
            std::cout << "Enter Left Recurse" << std::endl;
            std::vector<int> leftRecurs = AlignQueryToBWT_SuffixArray(ref, query, bwt, left, mid - 1, querySeqLength, querySeqIndex);

            std::cout << "Enter Right Recurse" << std::endl;
            std::vector<int> rightRecurs = AlignQueryToBWT_SuffixArray(ref, query,bwt, mid + 1, right, querySeqLength, querySeqIndex);

            // Add left and right recursive to matches
            for (int i = 0; i < leftRecurs.size(); i++) {
                matches.push_back(leftRecurs[i]);
            }
            for (int i = 0; i < rightRecurs.size(); i++) {
                matches.push_back(rightRecurs[i]);
            }

            return matches;
        }

        // recalculate midpoint
        mid = left + (right - left) / 2;
    }

    return matches;
}

BWTArray BWTFromSuffixArray(FReference &seq, int seqIndex) {
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

    // Determine number of unique characters
    //std::set<char> uniqueChars;
    for (int i = 0; i < bwt.last.length(); i++) {
        bwt.uniqueChars.insert(bwt.last[i]);
    }

    // Initialize counters
    bwt.counts.resize(bwt.uniqueChars.size());
    bwt.totalCounts.resize(bwt.uniqueChars.size());
    bwt.rankIndex.resize(bwt.uniqueChars.size());

    // Get total counts for each unique char
    int pos = 0;
    for (std::set<char>::iterator i = bwt.uniqueChars.begin(); i != bwt.uniqueChars.end(); i++) {
        bwt.totalCounts[pos] = std::count(bwt.last.begin(), bwt.last.end(), *i);
        pos++;
    }

    // Get rank index from total counts
    bwt.rankIndex[0] = bwt.last.length() - bwt.totalCounts[0];
    for (int i = bwt.uniqueChars.size() - 1; i > 0; i--) {
        if (i == bwt.uniqueChars.size() - 1) {
            bwt.rankIndex[i] = bwt.rankIndex[0] - bwt.totalCounts[i];
        }
        else {
            bwt.rankIndex[i] = bwt.rankIndex[i + 1] - bwt.totalCounts[i];
        }
    }

    for (int i = 0; i < bwt.last.length(); i++) {
        // find matching char, update counts and rank
        int pos = 0;
        for (std::set<char>::iterator j = bwt.uniqueChars.begin(); j != bwt.uniqueChars.end(); j++) {
            if (strncmp( bwt.last.substr(i, 1).c_str(), &*j, 1) == 0) {
                bwt.counts[pos]++;
                bwt.rank[i] = bwt.rankIndex[pos] + bwt.counts[pos];
                break;
            }
            pos++; // acts as j index
        }
    }
}

std::string FullSequenceFromBWTAndSuffixArray(FReference &seq, int seqIndex, BWTArray &bwt) {
    // Calculate rank
    std::cout << "Calculating Rank for BWT..." << std::endl;
    //CalculateBWTRank(bwt);
    std::cout << "Rank for BWT Calculated." << std::endl;

    // create output string 
    std::string output = "";
    output.reserve(bwt.last.length());
    int r = 0;
    
    // init
    output.append(bwt.last.substr(r,1));

    std::cout << "Recovering Sequence from BWT..." << std::endl;

    for (int i = 0; i < bwt.last.length(); i++) {
        //std::cout << "r: " << r << " c: " << bwt.last.substr(bwt.rank[r], 1).c_str() << std::endl; // Debug
        output.append(bwt.last.substr(bwt.rank[r], 1));
        r = bwt.rank[r];
    }
    // Fix missing character
    std::reverse(output.begin(), output.end());
    output.erase(0, 1);
    output.append("$");
    //std::cout << output << " " << output.length() << std::endl;

    std::cout << "Sequence Recovered." << std::endl;
    return output;
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

        std::cout << bwt.last << std::endl;

        // Calculate rank
        std::cout << "Calculating Rank for BWT..." << std::endl;
        CalculateBWTRank(bwt);
        std::cout << "Rank for BWT Calculated." << std::endl;

        // Recover original sequence from BWT and suffix Array and time it
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::cout << "Recovering Sequence from BWT and Suffix Array..." << std::endl;

        ref.Sequence.push_back({ 0 });
        ref.Sequence[0] = FullSequenceFromBWTAndSuffixArray(ref, 0, bwt);
        ref.SequenceLength.push_back(0);
        ref.SequenceLength[0] = ref.Sequence[0].length();
        ref.SequenceName.push_back("");
        ref.SequenceName[0] = "chr22";

        // Print and time output 
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto durationFPextraction = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "Sequence recovery time (ms): " << durationFPextraction << std::endl;

        // Load Queries
        std::cout << "Loading Queries..." << std::endl;
        FReference query(argv[3]);
        std::cout << "Queries Loaded." << std::endl;

        std::cout << "Aligning Queries..." << std::endl;
        
        
        std::cout << "=========== Sequence + Suffix Array ===========" << std::endl;
        std::vector<alignment> alignments;
        std::vector<int> alignmentIndexes;
        alignments.reserve(ref.Sequence.size());
        alignmentIndexes.reserve(ref.Sequence.size());

        std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
        // All queries
        for (int i = 0; i < query.Sequence.size(); i++) {
            // ref, query, left, right, querylength, query index, (ref index)
            alignmentIndexes = AlignQueryToSeq_SuffixArray(ref, query, 0, ref.SequenceLength[0] - 1, query.SequenceLength[i] - 1, i);

            // add new matches to alignment vector
            for (int j = 0; j < alignmentIndexes.size(); j++) {
                alignments.push_back(alignment(query.SequenceName[i], ref.SequenceName[0], alignmentIndexes[j]));
            }

            // Time step per iteration: end 
            //std::cout << "Alignment #" << i + 1 << std::endl;
            //std::cout << "Num Alignments Found: " << alignmentIndexes.size() << std::endl;
        }
        std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        auto durationFPextraction4 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
        std::cout << "Sequence + Suffix Array Alignment time Chr22 - 1,000,000 queries (ms): " << durationFPextraction4 << std::endl;

        std::cout << "Queries Aligned." << std::endl;

        std::cout << "=========== BWT + Suffix Array ===========" << std::endl;
        std::vector<alignment> alignments2;
        std::vector<int> alignmentIndexes2;
        alignments2.reserve(ref.Sequence.size());
        alignmentIndexes2.reserve(ref.Sequence.size());

        std::chrono::high_resolution_clock::time_point t7 = std::chrono::high_resolution_clock::now();
        // All queries
        for (int i = 0; i < query.Sequence.size(); i++) {
            // ref, query, left, right, querylength, query index, (ref index)
            alignmentIndexes2 = AlignQueryToBWT_SuffixArray(ref, query,bwt, 0, ref.SequenceLength[0] - 1, query.SequenceLength[i] - 1, i);

            // add new matches to alignment vector
            for (int j = 0; j < alignmentIndexes2.size(); j++) {
                alignments2.push_back(alignment(query.SequenceName[i], ref.SequenceName[0], alignmentIndexes2[j]));
            }

            // Time step per iteration: end 
            std::cout << "Alignment #" << i + 1 << std::endl;
            std::cout << "Num Alignments Found: " << alignmentIndexes2.size() << std::endl;
        }
        std::chrono::high_resolution_clock::time_point t8 = std::chrono::high_resolution_clock::now();
        auto durationFPextraction5 = std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count();
        std::cout << "BWT + Suffix Array Alignment time Chr22 - 1,000,000 queries (ms): " << durationFPextraction5 << std::endl;

        std::cout << "Queries Aligned (Again...)" << std::endl;


        // Print alignments to file
        //PrintAlignmentsToTerminal(alignments);

        // Save indexes to output file
        std::cout << "Saving alignments..." << std::endl;
        SaveAlignmentsToFile(GetFilename(argv[4]), alignments);
        SaveAlignmentsToFile("BWT_"+GetFilename(argv[4]), alignments2);

        std::cout << "Alignments saved." << std::endl;
        
        std::cout << "Program completed." << std::endl;
    }

    return 0;
}

