// Complete paired gzip FASTQ multi-length word/span evidence, C++17/zlib.
// Strict gzip/FASTQ parser adapted from the audited CODE8 scanner.
// Coordinates are zero-based half-open; pair_index is one-based.
// Read-pair support denotes distinct FASTQ records, not independent molecules.
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using U64 = std::uint64_t;



static void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
static bool identifier(const std::string& value) {
    if (value.empty() || value.size() > 80) return false;
    for (unsigned char ch : value)
        if (!((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ||
              (ch >= '0' && ch <= '9') || ch == '_' || ch == '.')) return false;
    return true;
}
static int base_index(char base) {
    switch (base) {case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1;}
}
static char complement(char base) {
    switch (base) {case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N';}
}
static std::string reverse_complement(const std::string& word) {
    std::string result(word.size(), 'N');
    for (std::size_t i = 0; i < word.size(); ++i) result[word.size()-i-1] = complement(word[i]);
    return result;
}
static std::vector<std::string> split(const std::string& line, char sep) {
    std::vector<std::string> answer;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= line.size(); ++i) if (i == line.size() || line[i] == sep) {
        answer.push_back(line.substr(start, i-start)); start = i+1;
    }
    return answer;
}
static U64 integer_argument(const std::string& text, U64 minimum, U64 maximum, const std::string& name) {
    require(!text.empty() && text.size() <= 20, "Invalid integer for " + name);
    U64 value = 0;
    for (char ch : text) {
        require(ch >= '0' && ch <= '9', "Invalid integer for " + name);
        require(value <= (maximum - U64(ch-'0')) / 10, "Integer exceeds bound for " + name);
        value = value*10 + U64(ch-'0');
    }
    require(value >= minimum && value <= maximum, "Integer outside bounds for " + name);
    return value;
}

// Stream gzip members directly with inflate, checking every CRC/trailer and
// rejecting non-gzip trailing bytes. Buffers are fixed at 64 KiB each.
class GzipLines {
    std::ifstream input;
    z_stream stream{};
    std::array<unsigned char, 65536> compressed{};
    std::array<unsigned char, 65536> decompressed{};
    std::size_t position = 0, available = 0;
    bool source_eof = false, member_ended = false, ended = false;
    std::string name;
    U64 members = 0;
    bool fill_input() {
        if (stream.avail_in) return true;
        if (source_eof) return false;
        input.read(reinterpret_cast<char*>(compressed.data()), compressed.size());
        auto count = input.gcount();
        require(!input.bad(), "Cannot read compressed input: " + name);
        if (count == 0) {source_eof = true; return false;}
        stream.next_in = compressed.data(); stream.avail_in = static_cast<uInt>(count);
        return true;
    }
    bool refill() {
        if (ended) return false;
        position = available = 0;
        while (available == 0) {
            if (member_ended) {
                if (!fill_input()) {ended = true; return false;}
                Bytef* next = stream.next_in; uInt remaining = stream.avail_in;
                require(inflateReset2(&stream, 15+16) == Z_OK, "Cannot reset gzip decoder");
                stream.next_in = next; stream.avail_in = remaining;
                member_ended = false;
            }
            require(fill_input(), "Truncated or missing gzip member: " + name);
            stream.next_out = decompressed.data(); stream.avail_out = decompressed.size();
            int status = inflate(&stream, Z_NO_FLUSH);
            available = decompressed.size() - stream.avail_out;
            if (status == Z_STREAM_END) {member_ended = true; ++members;}
            else if (status != Z_OK && status != Z_BUF_ERROR)
                throw std::runtime_error("Invalid gzip stream/CRC/trailing bytes: " + name + ": " +
                                         (stream.msg ? stream.msg : "inflate error"));
            if (available) return true;
            if (status == Z_BUF_ERROR && stream.avail_in)
                throw std::runtime_error("Gzip decoder made no progress: " + name);
        }
        return true;
    }
public:
    explicit GzipLines(const std::string& path): input(path, std::ios::binary), name(path) {
        require(bool(input), "Cannot open compressed input: " + path);
        require(inflateInit2(&stream, 15+16) == Z_OK, "Cannot initialize gzip decoder");
    }
    ~GzipLines() {inflateEnd(&stream);}
    GzipLines(const GzipLines&) = delete;
    GzipLines& operator=(const GzipLines&) = delete;
    bool line(std::string& result) {
        result.clear();
        while (true) {
            if (position == available && !refill()) {
                require(result.empty(), "FASTQ final line lacks newline: " + name);
                return false;
            }
            const auto begin = position;
            while (position < available && decompressed[position] != '\n') ++position;
            result.append(reinterpret_cast<const char*>(decompressed.data()+begin), position-begin);
            require(result.size() <= 4096, "FASTQ line exceeds 4096 bytes: " + name);
            if (position < available) {
                ++position;
                if (!result.empty() && result.back() == '\r') result.pop_back();
                return true;
            }
        }
    }
    bool complete() const {return ended;}
    U64 member_count() const {return members;}
};

struct Read {std::string id, sequence, quality;};
static std::string read_identifier(const std::string& text, int mate) {
    require(!text.empty(), "Empty FASTQ header identifier");
    for (unsigned char ch : text) if (!((ch >= 32 && ch <= 126) || ch == '\t'))
        throw std::runtime_error("Non-ASCII/invalid FASTQ header");
    std::istringstream tokens(text); std::string id, second; tokens >> id;
    require(!id.empty(), "Missing FASTQ identifier");
    if (id.size() >= 2 && id[id.size()-2] == '/' && (id.back() == '1' || id.back() == '2')) {
        require(id.back()-'0' == mate, "FASTQ mate tag differs from input mate"); id.resize(id.size()-2);
    }
    require(!id.empty(), "Missing normalized FASTQ identifier");
    if (tokens >> second) {
        std::size_t digits = 0;
        while (digits < second.size() && second[digits] >= '0' && second[digits] <= '9') ++digits;
        if (digits && digits < second.size() && second[digits] == ':') {
            auto second_colon=second.find(':',digits+1);
            auto third_colon=second_colon==std::string::npos ? std::string::npos : second.find(':',second_colon+1);
            bool valid=digits==1 && (second[0]=='1' || second[0]=='2') && second_colon==3 &&
                       (second[2]=='Y' || second[2]=='N') && third_colon!=std::string::npos &&
                       third_colon>second_colon+1 && third_colon+1<second.size();
            if(valid) for(std::size_t i=second_colon+1;i<third_colon;++i)
                if(second[i]<'0' || second[i]>'9') valid=false;
            require(valid, "Malformed CASAVA mate header");
            require(second[0]-'0' == mate, "CASAVA mate tag differs from input mate");
        }
    }
    return id;
}
static bool next_read(GzipLines& file, int mate, Read& read, U64 pair_index) {
    std::string header, separator;
    if (!file.line(header)) return false;
    const std::string where = " at pair " + std::to_string(pair_index) + " mate " + std::to_string(mate);
    require(!header.empty() && header[0] == '@', "FASTQ header must begin with @" + where);
    read.id = read_identifier(header.substr(1), mate);
    require(file.line(read.sequence) && file.line(separator) && file.line(read.quality), "Truncated FASTQ record" + where);
    require(!separator.empty() && separator[0] == '+', "FASTQ separator must begin with +" + where);
    if (separator.size() > 1) require(read_identifier(separator.substr(1), mate) == read.id, "Repeated + identifier differs" + where);
    require(!read.sequence.empty() && read.sequence.size() <= 512, "Read length must be 1..512" + where);
    for (char& ch : read.sequence) {
        if (ch >= 'a' && ch <= 'z') ch = char(ch-'a'+'A');
        switch(ch) {
            case 'A': case 'C': case 'G': case 'T': case 'R': case 'Y': case 'S':
            case 'W': case 'K': case 'M': case 'B': case 'D': case 'H': case 'V': case 'N': break;
            default: throw std::runtime_error("Non-IUPAC sequence base" + where);
        }
    }
    require(read.sequence.size() == read.quality.size(), "Sequence/quality lengths differ" + where);
    for (unsigned char q : read.quality) if(q < 33 || q > 126)
        throw std::runtime_error("Invalid Phred+33 quality" + where);
    return true;
}

static std::vector<std::vector<std::string>> table(const fs::path& path, const std::vector<std::string>& header) {
    require(fs::is_regular_file(path) && !fs::is_symlink(path), "Catalog must be a regular non-symlink file");
    require(fs::file_size(path) <= 64*1024*1024, "Catalog exceeds 64 MiB bound");
    std::ifstream file(path); std::string line;
    require(bool(std::getline(file, line)), "Empty catalog: " + path.string());
    if (!line.empty() && line.back() == '\r') line.pop_back();
    require(split(line, '\t') == header, "Unexpected TSV header: " + path.string());
    std::vector<std::vector<std::string>> rows;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        require(!line.empty() && line.size() <= 4096, "Empty/oversized catalog row");
        auto row = split(line, '\t'); require(row.size() == header.size(), "Wrong number of catalog columns");
        rows.push_back(std::move(row)); require(rows.size() <= 150000, "Too many catalog rows");
    }
    require(!file.bad(), "Cannot read catalog"); return rows;
}
struct Query {std::string id,sequence;};
struct Group {std::string id,state;int k,span;std::vector<int> words;};
struct Hit {
    int start,query;bool q20;
    bool operator<(const Hit& b) const {return std::tie(start,query,q20)<std::tie(b.start,b.query,b.q20);}
    bool operator==(const Hit& b) const {return start==b.start&&query==b.query&&q20==b.q20;}
};
struct Node {std::array<int,4> next{{-1,-1,-1,-1}};int fail=0;std::vector<int> output;};
class Automaton {
    std::vector<Node> nodes{1};
    const std::vector<Query>& queries;
public:
    explicit Automaton(const std::vector<Query>& q):queries(q) {
        for(std::size_t i=0;i<q.size();++i) {
            auto rev=reverse_complement(q[i].sequence);
            for(int strand=0;strand<2;++strand) {
                if(strand==1&&rev==q[i].sequence) continue;
                const auto& seq=strand?rev:q[i].sequence;int state=0;
                for(char c:seq) {
                    int b=base_index(c),next=nodes[state].next[b];
                    if(next<0) {next=static_cast<int>(nodes.size());nodes[state].next[b]=next;nodes.emplace_back();}
                    state=next;
                }
                nodes[state].output.push_back(static_cast<int>(i));
            }
        }
        require(nodes.size()<=10000000,"Automaton state cap exceeded");
        std::queue<int> queue;
        for(int b=0;b<4;++b) {
            int child=nodes[0].next[b];
            if(child<0) nodes[0].next[b]=0;else queue.push(child);
        }
        while(!queue.empty()) {
            int state=queue.front();queue.pop();
            for(int b=0;b<4;++b) {
                int child=nodes[state].next[b];
                if(child<0) {nodes[state].next[b]=nodes[nodes[state].fail].next[b];continue;}
                int failure=nodes[nodes[state].fail].next[b];nodes[child].fail=failure;
                nodes[child].output.insert(nodes[child].output.end(),nodes[failure].output.begin(),nodes[failure].output.end());
                queue.push(child);
            }
        }
    }
    std::vector<Hit> scan(const Read& read) const {
        std::vector<Hit> hits;int state=0,last_low=-1;
        for(std::size_t i=0;i<read.sequence.size();++i) {
            if(static_cast<unsigned char>(read.quality[i])<53) last_low=static_cast<int>(i);
            int b=base_index(read.sequence[i]);
            if(b<0) {state=0;continue;}
            state=nodes[state].next[b];
            for(int query:nodes[state].output) {
                int start=static_cast<int>(i)+1-static_cast<int>(queries[query].sequence.size());
                hits.push_back({start,query,last_low<start});
            }
        }
        return hits;
    }
    std::size_t size() const {return nodes.size();}
};
struct DirectPattern {std::string sequence;int query;};
static std::vector<DirectPattern> direct_patterns(const std::vector<Query>& queries) {
    std::vector<DirectPattern> patterns;
    patterns.reserve(2*queries.size());
    for(std::size_t i=0;i<queries.size();++i) {
        patterns.push_back({queries[i].sequence,static_cast<int>(i)});
        auto reverse=reverse_complement(queries[i].sequence);
        if(reverse!=queries[i].sequence) patterns.push_back({std::move(reverse),static_cast<int>(i)});
    }
    return patterns;
}
static std::vector<Hit> direct_scan(const Read& read,const std::vector<DirectPattern>& patterns) {
    std::vector<Hit> hits;
    for(const auto& pattern:patterns) {
        const auto& word=pattern.sequence;
        auto start=read.sequence.find(word);
        while(start!=std::string::npos) {
            bool good=true;
            for(std::size_t j=start;j<start+word.size();++j) if(static_cast<unsigned char>(read.quality[j])<53) good=false;
            hits.push_back({static_cast<int>(start),pattern.query,good});start=read.sequence.find(word,start+1);
        }
    }
    std::sort(hits.begin(),hits.end());return hits;
}
struct WordStats {U64 raw_occurrences=0,q20_occurrences=0,raw_read_support=0,q20_read_support=0,raw_pair_support=0,q20_pair_support=0;};
struct GroupStats {U64 raw_same_read_complete_reads=0,q20_same_read_complete_reads=0,raw_same_read_complete_pairs=0,q20_same_read_complete_pairs=0,raw_mate_union_complete_pairs=0,q20_mate_union_complete_pairs=0;};
class Counters {
    const std::vector<Group>& groups;
    const std::vector<std::vector<int>>& owners;
    std::vector<U64> raw_read_seen,q20_read_seen,raw_pair_seen,q20_pair_seen;
    std::vector<U64> raw_read_group_seen,q20_read_group_seen,raw_pair_group_seen,q20_pair_group_seen;
    std::vector<U64> raw_read_group_count,q20_read_group_count,raw_pair_group_count,q20_pair_group_count;
    std::vector<U64> raw_same_pair_seen,q20_same_pair_seen;
public:
    std::vector<WordStats> words;
    std::vector<GroupStats> group_counts;
    std::vector<int> touched;
    Counters(std::size_t n,const std::vector<Group>& g,const std::vector<std::vector<int>>& o):
        groups(g),owners(o),raw_read_seen(n),q20_read_seen(n),raw_pair_seen(n),q20_pair_seen(n),
        raw_read_group_seen(g.size()),q20_read_group_seen(g.size()),raw_pair_group_seen(g.size()),q20_pair_group_seen(g.size()),
        raw_read_group_count(g.size()),q20_read_group_count(g.size()),raw_pair_group_count(g.size()),q20_pair_group_count(g.size()),
        raw_same_pair_seen(g.size()),q20_same_pair_seen(g.size()),words(n),group_counts(g.size()) {}
    void add_read(const std::vector<Hit>& hits,U64 pair,int mate) {
        U64 read=pair*2+mate;
        for(const auto& h:hits) {
            int q=h.query;auto& w=words[q];++w.raw_occurrences;if(h.q20) ++w.q20_occurrences;
            bool new_raw_read=raw_read_seen[q]!=read,new_raw_pair=raw_pair_seen[q]!=pair;
            bool new_q20_read=h.q20&&q20_read_seen[q]!=read,new_q20_pair=h.q20&&q20_pair_seen[q]!=pair;
            if(new_raw_read) {raw_read_seen[q]=read;++w.raw_read_support;}
            if(new_raw_pair) {raw_pair_seen[q]=pair;++w.raw_pair_support;}
            if(new_q20_read) {q20_read_seen[q]=read;++w.q20_read_support;}
            if(new_q20_pair) {q20_pair_seen[q]=pair;++w.q20_pair_support;}
            if(!(new_raw_read||new_raw_pair||new_q20_read||new_q20_pair)) continue;
            for(int g:owners[q]) {
                auto required=groups[g].words.size();auto& c=group_counts[g];
                if(new_raw_read) {
                    if(raw_read_group_seen[g]!=read) {raw_read_group_seen[g]=read;raw_read_group_count[g]=0;}
                    if(++raw_read_group_count[g]==required) {
                        ++c.raw_same_read_complete_reads;
                        if(raw_same_pair_seen[g]!=pair) {raw_same_pair_seen[g]=pair;++c.raw_same_read_complete_pairs;}
                    }
                }
                if(new_raw_pair) {
                    if(raw_pair_group_seen[g]!=pair) {raw_pair_group_seen[g]=pair;raw_pair_group_count[g]=0;touched.push_back(g);}
                    if(++raw_pair_group_count[g]==required) ++c.raw_mate_union_complete_pairs;
                }
                if(new_q20_read) {
                    if(q20_read_group_seen[g]!=read) {q20_read_group_seen[g]=read;q20_read_group_count[g]=0;}
                    if(++q20_read_group_count[g]==required) {
                        ++c.q20_same_read_complete_reads;
                        if(q20_same_pair_seen[g]!=pair) {q20_same_pair_seen[g]=pair;++c.q20_same_read_complete_pairs;}
                    }
                }
                if(new_q20_pair) {
                    if(q20_pair_group_seen[g]!=pair) {q20_pair_group_seen[g]=pair;q20_pair_group_count[g]=0;}
                    if(++q20_pair_group_count[g]==required) ++c.q20_mate_union_complete_pairs;
                }
            }
        }
    }
    void add_pair(const std::vector<Hit>& a,const std::vector<Hit>& b,U64 pair) {
        touched.clear();add_read(a,pair,0);add_read(b,pair,1);
    }
    // Four mutually exclusive evidence categories per quality mode.
    int category(int g,U64 pair,bool quality) const {
        const auto& seen=quality?q20_pair_group_seen:raw_pair_group_seen;
        if(seen[g]!=pair) return -1;
        const auto& fullseen=quality?q20_pair_seen:raw_pair_seen;
        if(fullseen[groups[g].span]==pair) return 2;
        const auto& sameseen=quality?q20_same_pair_seen:raw_same_pair_seen;
        if(sameseen[g]==pair) return 0;
        const auto& counts=quality?q20_pair_group_count:raw_pair_group_count;
        return counts[g]==groups[g].words.size()?1:3;
    }
};
static const char* RECORD_HEADER="pair_index\tread_id\tr1_sequence\tr1_quality\tr2_sequence\tr2_quality\n";
static void record(std::ofstream& out,U64 index,const Read& a,const Read& b) {
    out<<index<<'\t'<<a.id<<'\t'<<a.sequence<<'\t'<<a.quality<<'\t'<<b.sequence<<'\t'<<b.quality<<'\n';
    require(bool(out),"Cannot write bounded FASTQ examples");
}
static void write_counts(const fs::path& path,const Counters& counts,const std::vector<Query>& queries,const std::vector<Group>& groups,const std::string& prefix) {
    std::ofstream w(path/(prefix+"word_counts.tsv.partial"));
    w<<"pattern_id\tsequence\traw_occurrences\tq20_occurrences\traw_read_support\tq20_read_support\traw_pair_support\tq20_pair_support\n";
    for(std::size_t i=0;i<queries.size();++i) {
        const auto& x=counts.words[i];w<<queries[i].id<<'\t'<<queries[i].sequence<<'\t'<<x.raw_occurrences<<'\t'<<x.q20_occurrences<<'\t'<<x.raw_read_support<<'\t'<<x.q20_read_support<<'\t'<<x.raw_pair_support<<'\t'<<x.q20_pair_support<<'\n';
    }
    w.close();require(bool(w),"Word count output failed");
    std::ofstream g(path/(prefix+"group_counts.tsv.partial"));
    g<<"group_id\tstate_id\tk\traw_same_read_complete_reads\tq20_same_read_complete_reads\traw_same_read_complete_pairs\tq20_same_read_complete_pairs\traw_mate_union_complete_pairs\tq20_mate_union_complete_pairs\traw_exact_span_pairs\tq20_exact_span_pairs\n";
    for(std::size_t i=0;i<groups.size();++i) {
        const auto& x=counts.group_counts[i];const auto& q=counts.words[groups[i].span];
        g<<groups[i].id<<'\t'<<groups[i].state<<'\t'<<groups[i].k<<'\t'<<x.raw_same_read_complete_reads<<'\t'<<x.q20_same_read_complete_reads<<'\t'<<x.raw_same_read_complete_pairs<<'\t'<<x.q20_same_read_complete_pairs<<'\t'<<x.raw_mate_union_complete_pairs<<'\t'<<x.q20_mate_union_complete_pairs<<'\t'<<q.raw_pair_support<<'\t'<<q.q20_pair_support<<'\n';
    }
    g.close();require(bool(g),"Group count output failed");
}
int main(int argc,char** argv) {
    try {
        std::map<std::string,std::string> args;
        for(int i=1;i<argc;i+=2) {
            require(i+1<argc,"All options require values");std::string key=argv[i];
            require(std::set<std::string>{"--queries","--groups","--r1","--r2","--out","--validation-pairs","--examples-per-category"}.count(key),"Unknown option");
            require(args.emplace(key,argv[i+1]).second,"Duplicate option");
        }
        for(const auto& key:{"--queries","--groups","--r1","--r2","--out"}) require(args.count(key),std::string("Missing ")+key);
        U64 validation=integer_argument(args.count("--validation-pairs")?args["--validation-pairs"]:"1000",1,10000,"validation-pairs");
        U64 examples=integer_argument(args.count("--examples-per-category")?args["--examples-per-category"]:"3",1,10,"examples-per-category");
        std::vector<Query> queries;std::map<std::string,int> ids;std::set<std::string> sequences;
        for(const auto& row:table(args["--queries"],{"pattern_id","sequence"})) {
            require(identifier(row[0])&&ids.emplace(row[0],static_cast<int>(queries.size())).second,"Invalid/duplicate query ID");
            require(row[1].size()>=1&&row[1].size()<=512,"Query length must be1..512");
            for(char c:row[1]) require(base_index(c)>=0,"Queries require uppercase ACGT");
            require(row[1]<=reverse_complement(row[1])&&sequences.insert(row[1]).second,"Queries must be unique and canonical");
            queries.push_back({row[0],row[1]});
        }
        require(!queries.empty()&&queries.size()<=150000,"Query count must be1..150000");
        std::vector<Group> groups;std::set<std::string> gids;std::set<std::pair<std::string,int>> groupkeys;
        std::vector<std::vector<int>> owners(queries.size());
        for(const auto& row:table(args["--groups"],{"group_id","state_id","k","word_ids","span_query_id"})) {
            require(identifier(row[0])&&gids.insert(row[0]).second&&identifier(row[1]),"Invalid/duplicate group identity");
            int k=static_cast<int>(integer_argument(row[2],1,512,"group k"));
            require(groupkeys.emplace(row[1],k).second&&ids.count(row[4]),"Duplicate group or unknown span query");
            int span=ids.at(row[4]);require(k<=static_cast<int>(queries[span].sequence.size()),"k exceeds span length");
            Group group{row[0],row[1],k,span,{}};std::set<std::string> found,expected;
            for(const auto& id:split(row[3],',')) {
                require(ids.count(id),"Unknown group word query");int q=ids.at(id);
                require(static_cast<int>(queries[q].sequence.size())==k&&found.insert(queries[q].sequence).second,"Wrong length or duplicate group word");
                group.words.push_back(q);owners[q].push_back(static_cast<int>(groups.size()));
            }
            const auto& seq=queries[span].sequence;
            for(std::size_t i=0;i+k<=seq.size();++i) {auto word=seq.substr(i,k);expected.insert(std::min(word,reverse_complement(word)));}
            require(found==expected,"Group word queries must equal every canonical span window");groups.push_back(std::move(group));
        }
        require(!groups.empty()&&groups.size()<=5000,"Group count must be1..5000");
        Automaton automaton(queries);auto direct=direct_patterns(queries);Counters counts(queries.size(),groups,owners),prefix_counts(queries.size(),groups,owners);
        fs::path out=args["--out"];require(!fs::exists(out)&&fs::create_directory(out),"Output must be a new directory");
        std::ofstream prefix(out/"validation_pairs.tsv.partial",std::ios::binary),records(out/"example_records.tsv.partial",std::ios::binary),selection(out/"example_selections.tsv.partial");
        require(bool(prefix)&&bool(records)&&bool(selection),"Cannot open native evidence output");
        prefix<<RECORD_HEADER;records<<RECORD_HEADER;selection<<"group_id\tcategory\tpair_index\n";
        std::vector<std::array<U64,8>> selected(groups.size());
        const std::array<std::string,4> category_names={"same_read_complete_no_exact_span","mate_union_only","exact_span","partial_only"};
        GzipLines r1(args["--r1"]),r2(args["--r2"]);Read a,b;U64 pairs=0,bases=0,samples=0,selections=0;double oracle_seconds=0;auto start=Clock::now();
        while(true) {
            bool first=next_read(r1,1,a,pairs+1),second=next_read(r2,2,b,pairs+1);
            require(first==second,"FASTQ mates have unequal record counts");if(!first) break;
            require(a.id==b.id,"FASTQ paired read identifiers disagree");++pairs;
            require(pairs<=1000000000000ULL,"Input pair count exceeds uint64 safety bound");bases+=a.sequence.size()+b.sequence.size();
            auto x=automaton.scan(a),y=automaton.scan(b);counts.add_pair(x,y,pairs);
            if(pairs<=validation) {
                auto oracle_start=Clock::now();auto sx=x,sy=y;std::sort(sx.begin(),sx.end());std::sort(sy.begin(),sy.end());
                require(sx==direct_scan(a,direct)&&sy==direct_scan(b,direct),"Independent direct-find oracle mismatch");
                oracle_seconds+=std::chrono::duration<double>(Clock::now()-oracle_start).count();
                prefix_counts.add_pair(x,y,pairs);record(prefix,pairs,a,b);
            }
            bool retain=false;
            for(int g:counts.touched) for(int q=0;q<2;++q) {
                int category=counts.category(g,pairs,q!=0);if(category<0) continue;int slot=q*4+category;
                if(selected[g][slot]>=examples) continue;
                ++selected[g][slot];++selections;retain=true;
                selection<<groups[g].id<<'\t'<<(q?"q20_":"raw_")<<category_names[category]<<'\t'<<pairs<<'\n';
            }
            if(retain) {record(records,pairs,a,b);++samples;}
            if(pairs%1000000==0) std::cerr<<"PROGRESS pairs="<<pairs<<" examples="<<samples<<"\n";
        }
        require(pairs>0&&r1.complete()&&r2.complete(),"No paired records or incomplete gzip EOF");
        prefix.close();records.close();selection.close();require(bool(prefix)&&bool(records)&&bool(selection),"Evidence output write failed");
        write_counts(out,counts,queries,groups,"");write_counts(out,prefix_counts,queries,groups,"prefix_");
        std::ofstream summary(out/"summary.json.partial");
        summary<<std::setprecision(12)<<"{\"status\":\"PANEL_SCAN_COMPLETE\",\"schema_version\":1,\"full_eof\":true,\"min_phred\":20,\"pairs\":"<<pairs<<",\"reads\":"<<2*pairs<<",\"bases\":"<<bases<<",\"queries\":"<<queries.size()<<",\"groups\":"<<groups.size()<<",\"validation_pairs\":"<<std::min(pairs,validation)<<",\"requested_validation_pairs\":"<<validation<<",\"oracle_mismatches\":0,\"example_records\":"<<samples<<",\"example_selections\":"<<selections<<",\"examples_per_category\":"<<examples<<",\"r1_gzip_members\":"<<r1.member_count()<<",\"r2_gzip_members\":"<<r2.member_count()<<",\"automaton_states\":"<<automaton.size()<<",\"direct_oracle_seconds\":"<<oracle_seconds<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<"}\n";
        summary.close();require(bool(summary),"Summary write failed");
        for(const auto& name:{"validation_pairs.tsv","example_records.tsv","example_selections.tsv","word_counts.tsv","group_counts.tsv","prefix_word_counts.tsv","prefix_group_counts.tsv","summary.json"}) fs::rename(out/(std::string(name)+".partial"),out/name);
        std::cout<<"PANEL_SCAN_COMPLETE pairs="<<pairs<<" example_records="<<samples<<"\n";return 0;
    } catch(const std::exception& e) {std::cerr<<"PANEL_SCAN_ERROR "<<e.what()<<"\n";return 1;}
}
