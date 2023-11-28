/*
 * chi2bgen.cc
 */

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <zlib.h>
#include <stdlib.h>
#include <omp.h>
#include <boost/log/common.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>
#include <boost/bind/bind.hpp>
#include <iostream>
#include <chrono>
#include <ctime>

#include "genfile/bgen/bgen.hpp"
#include "genfile/types.hpp"

#include "chi/chiReader.hpp"
// #include "bazaar/Index.hh"

#include "version.h"
#include <iplogger.h>

namespace logging = boost::log;

std::string usage =
    "\n"
    "Converts a .chi file to a bgen file\n"
    "\n"
    "USAGE:\n"
    "\tchi2bgen [-p | -u] [-n] [-b <prob_bits>] [-t <num_threads>] [-c <chunk_size> -i <index>] -f <chi_fileroot> -o <bgen_outfile>\n"
    "\tchi2bgen -l <bgen_infile>\n"
    "\tchi2bgen -h | --help\n"
    "\tchi2bgen -v | --version\n"
    "\n"
    "DESCRIPTION:\n"
    "\tchi2bgen converts the specified chi files into a bgen file.\n"
    "\n"
    "OPTIONS:\n"
    "\t-b=<prob_bits>\t\tHow many bits per genotype probability. Defaults to 8.\n"
    "\t-t=<num_threads>\tHow many threads to use for execution. Defaults to 4.\n"
    "\t-m=<mb>\t\t\tHow much RAM in GB to utilize for buffering. Defaults to 2.\n"
    "\t-c=<chunk_size>\t\tMarkers per chunk to process. Defaults to whole file.\n"
    "\t-i=<index>\t\tIndex of the chunk to process.\n"
    "\t-f=<chi_fileroot>\tPath to set of chi files (.pns, .markers, .chi) to convert.\n"
    "\t-o=<out_fileroot>\tPath to outfiles. (<out_fileroot>.bgen, <out_fileroot>.sample). Will be created.\n"
    "\t-l=<bgen_infile>\tRead provided bgenfile to standard out. Mostly for debugging purposes.\n"
    "\t-d --debug\t\tRun verbosely.\n"
    "\t-h --help\t\tPrint help and exit.\n"
    "\t-v --version\t\tPrint version and exit.\n"
    "\n"
    "EXAMPLES:\n"
    "\tchi2bgen -f /nfs/prog/bioinfo/extdata/deCODE/testinput/20170519_biotools/format_converters/test -o <PathToOutfiles>. \n"
    "\n";



// Dead simple option parsing
// -----------------------------
std::string getCmdOption(char** begin, char** end, const std::string& option) {

    char **itr = std::find(begin, end, option);

    if (itr != end && ++itr != end) {
        return std::string(*itr);
    }

    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {

    return std::find(begin, end, option) != end;

}
// -----------------------------


// Callback functions for bgen c++ implementation
// -----------------------------
std::string chi_get_tbx_allele(TabixIndex::Iterator& m, std::size_t i) {
    if (i == 0) {
        std::string A0_allele(m->field(3));
        return A0_allele;
    }
    else {
        std::string A1_allele(m->field(4));
        return A1_allele;
    }
}

std::string chi_get_allele(ImputedMarker* m, std::size_t i ) {

    if (i == 0) {
        std::string A0_allele(m->allele);
        return A0_allele;
    }
    else {
        std::string A1_allele(m->other_allele);
        return A1_allele;
    }
}


void samplevector_addid(std::vector<std::string>* sids, std::string id) {
    sids->push_back(id);
}
// -----------------------------


// Standard out setter - Needed to work with parse_probability_data API
// -----------------------------

struct StandardOutSetter {

    void initialise(std::size_t num_samples, std::size_t num_alleles) {

        std::cout << "Sample count: " << num_samples << "\n";
        std::cout << "Allele count: " << num_alleles << "\n";
    }

    void set_min_max_ploidy(uint32_t, uint32_t, uint32_t, uint32_t) {

        std::cout << "Redundant data :)\n";
    }

    bool set_sample(std::size_t i) {
        std::cout << "Reporting sample " << i << "\n";
        return true;
    }

    void set_number_of_entries(
        uint32_t ploidy_at_sample,
        uint32_t num_prob_values,
        genfile::OrderType ordert,
        genfile::ValueType valuet) {

            std::cout << "Ploidy at sample: " << ploidy_at_sample << "\n";
            std::cout << "Prob values for sample " << num_prob_values << "\n";
            std::cout << "OrderType is " << ordert << "\n";
            std::cout << "ValueType is " << valuet << "\n";
        }

    void set_value(std::size_t i, genfile::MissingValue) {
        std::cout << i << ": Missing data.";
    }

    void set_value(std::size_t i, double v) {
        std::cout << i << ": " << v << "\n";
    }

    void finalise() {
        std::cout << "--------\n";
    }
};
// -----------------------------


// Reads the given file and reports data to standard out.
// -----------------------------
void test_read(std::string testfile) {

    std::ifstream verify;
    verify.open(testfile);

    genfile::bgen::Context generated_context = genfile::bgen::Context();

    uint32_t offset;

    genfile::bgen::read_offset(verify, &offset);

    std::cout << "first offset: " << offset << ".\n";

    std::cout << "header size: " << genfile::bgen::read_header_block(verify, &generated_context) << ".\n";

    std::size_t written_marker_count = generated_context.number_of_variants;

    std::cout << "PN count : " << generated_context.number_of_samples << ", marker count: " << generated_context.number_of_variants << ".\n";
    std::cout << "flags: " << generated_context.flags << ".\n";

    // read sample identification
    std::vector<std::string> sample_ids;

    if (generated_context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(verify, generated_context, boost::bind(&samplevector_addid, &sample_ids, boost::placeholders::_1));
    }

    for (auto sid: sample_ids) {
        std::cout << sid << "\n";
    }

    std::string written_snpid;
    std::string written_rsid;
    std::string written_chromosome;
    uint32_t written_pos;
    std::string written_first_allele;
    std::string written_second_allele;

    std::vector<genfile::byte_t> buffer1, buffer2;

    StandardOutSetter myboy;

    for (size_t i = 0; i < written_marker_count; i++) {
        genfile::bgen::read_snp_identifying_data(
            verify,
            generated_context,
            &written_snpid,
            &written_rsid,
            &written_chromosome,
            &written_pos,
            &written_first_allele,
            &written_second_allele
        );

        std::cout << "Marker: " << i << " : Chromosome: " << written_chromosome << " Position: " << written_pos << "\n";
        std::cout << "Name: " << written_snpid << " RSID: " << written_rsid << "\n";
        std::cout << "First allele: " << written_first_allele << " Second allele: " << written_second_allele << "\n";

         genfile::bgen::read_and_parse_genotype_data_block(verify, generated_context, myboy, &buffer1, &buffer2);
    }
}
// -----------------------------
// Sample files are useful for downstream work and needed by plink and friends.

void create_sample_file(std::string outroot, std::vector<std::string> sample_ids) {
    std::ofstream sample_file;

    sample_file.open(outroot + ".sample");

    // Sample files contain a header line explaining the columns
    // Minimally acceptable is:
    sample_file << "ID_1 ID_2 missing\n";

    // sample files also have a second line
    // minimally acceptable is:
    sample_file << "0 0 0\n";

    for(auto sample_id: sample_ids) {
        sample_file << sample_id << " " << sample_id << " 0\n";
    }

    std::cout << "Wrote sample file " << outroot << ".sample\n";
}


// -----------------------------
// Sometimes users would prefer not having phasing information.
std::tuple<double, double, double> unphase(double pat_alt, double mat_alt) {
    // AA prob
    double ref_hom_prob = (1.0 - pat_alt) * (1.0 - mat_alt);
    // Aa/aA prob
    double het_prob = (((1.0 - pat_alt) * mat_alt) + (pat_alt * (1.0 - mat_alt)));
    // aa prob
    double alt_hom_prob = (pat_alt * mat_alt);

    return std::make_tuple(ref_hom_prob, het_prob, alt_hom_prob);
}


// -----------------------------
// Lets look for index files
std::string get_indexfname(const std::string &gt_filepath) {
    std::string indexfname;
    const std::string cmap1_base = gt_filepath + ".cmap1";
    const std::string cmap_base = gt_filepath + ".cmap";

    if (boost::filesystem::exists(cmap1_base + ".irng.gz")) {
        indexfname = cmap1_base + ".irng.gz";
    }
    else if (boost::filesystem::exists(cmap_base + ".gz")) {
        indexfname = cmap_base + ".irng.gz";
    } 
    else {
        BOOST_LOG_TRIVIAL(debug) << "Neither cmap1 nor cmap file found. Exiting. \
        Are the index files missing?" << std::endl;
        exit(1);
    }

    return indexfname;
}


// -------------------
// Initialize context.
genfile::bgen::Context init_context(
    int num_samples,
    int num_variants) {
    
    genfile::bgen::Context mycontext = genfile::bgen::Context();
    mycontext.number_of_samples = num_samples;
    mycontext.number_of_variants = num_variants;
    // Always sample id blocks.
    mycontext.flags |= genfile::bgen::e_SampleIdentifiers;
    // Always layout 2, zstd compression
    mycontext.flags |= genfile::bgen::e_Layout2;
    mycontext.flags = (mycontext.flags & ~genfile::bgen::e_CompressedSNPBlocks) | genfile::bgen::e_ZstdCompression;

    return mycontext;
}

// ------------------------------------------
// Write BGEN file header with given context.
// Overload for now, overrides old func probs.
void write_bgen_header(
    std::ofstream &bgenfile,
    genfile::bgen::Context context,
    std::vector<std::string> &pns) {
    
    uint32_t offset = context.header_size();
    genfile::bgen::write_offset(bgenfile, offset);
    // Will track offset, repeatedly writing, overwriting.
    // Thus it is done in BGEN example implementation, thus do we do.
    genfile::bgen::write_header_block(bgenfile, context);
    offset += genfile::bgen::write_sample_identifier_block(bgenfile, context, pns);
    bgenfile.seekp(0, std::ios_base::beg);
    genfile::bgen::write_offset(bgenfile, offset);
    genfile::bgen::write_header_block(bgenfile, context);
    // Seek bgenfile to start of variant data.
    bgenfile.seekp(offset + 4, std::ios_base::beg);
}

// -----------------------------
// Write bgen header for threads.
// Expect ofstream to be open.
// Not yet sure if unsigned or signed, leave as is until compiler complains.
genfile::bgen::Context write_bgen_header(
    std::ofstream &bgenfile, 
    int num_samples,
    int num_variants, 
    std::vector<std::string> &pns) {
    // The context tracks some values for the bgenfile header.
    genfile::bgen::Context mycontext = genfile::bgen::Context();
    mycontext.number_of_samples = num_samples;
    mycontext.number_of_variants = num_variants;

    uint32_t myoffset = mycontext.header_size();
    // Bgen code uses this pattern of repeatedly rewriting offset,
    // dependent on data written.
    // Which feels really silly. 
    // I'll avoid becoming an expert and just cargo-cult.
    genfile::bgen::write_offset(bgenfile, myoffset);
    // We will always include sample identifier blocks
    mycontext.flags |= genfile::bgen::e_SampleIdentifiers;
    // We will always use layout 2 with zstd compression.
    // Remember to remove compress option
    mycontext.flags |= genfile::bgen::e_Layout2;
    mycontext.flags = (mycontext.flags & ~genfile::bgen::e_CompressedSNPBlocks) | genfile::bgen::e_ZstdCompression;
    
    // With context set, we write it into the header block.
    genfile::bgen::write_header_block(bgenfile, mycontext);

    // Then we do the sample id block.
    myoffset += genfile::bgen::write_sample_identifier_block(bgenfile, mycontext, pns);
    // Must rewrite header offset after addition of sample id block
    bgenfile.seekp(0, std::ios_base::beg);
    genfile::bgen::write_offset(bgenfile, myoffset);
    genfile::bgen::write_header_block(bgenfile, mycontext);

    // With header block complete, we move the filepointer so caller
    // may just start writing variant data.
    bgenfile.seekp(myoffset + 4, std::ios_base::beg);

    return mycontext;
}



// -----------------------------
// Fill a writer object for a block of variant data,
void encode_block(
    float* vals, 
    size_t pn_count,
    genfile::bgen::GenotypeDataBlockWriter* writer) {
    
    genfile::MissingValue missing_allele;
    uint16_t allele_count = 2;

    writer->initialise(pn_count, allele_count);

    for (std::size_t i = 0; i < pn_count*2; i += 2) {
        writer->set_sample(i/2);
        writer->set_number_of_entries(2, 4, genfile::ePerPhasedHaplotypePerAllele, genfile::eProbability);

            // Every (n, n+1) where n=2k, pair in phased BGEN format data must sum to 1

            // data is present
            if (vals[i] != -1.0 && vals[i+1] != -1.0 ) {

                // paternal is reference
                writer->set_value(0, (1 - vals[i]));
                // paternal is alternate
                writer->set_value(1, vals[i]);

                // maternal is reference
                writer->set_value(2, (1 - vals[i+1]));
                // maternal is alternate
                writer->set_value(3, vals[i+1]);
            }
            // data is missing
            else {
                writer->set_value(0, missing_allele);
                writer->set_value(1, missing_allele);
                writer->set_value(2, missing_allele);
                writer->set_value(3, missing_allele);
            }            
    }
    writer->finalise();
} 


// -----------------------------


int main(int argc, char* argv[]) {

    ip_logger(my_version, argv, argc);
    std::ios_base::sync_with_stdio();
    /* ---------------------------- */
    // Command line parsing begins.
    /* ---------------------------- */
    // Parse everything from command line. Use a lib if program becomes more complicated.
    if (argc == 1) {
	std::cout << usage;
        return 0;
    }

    if (cmdOptionExists(argv, argv + argc, "-h") || cmdOptionExists(argv, argv + argc, "--help")) {
        // Write help string and exit
        std::cout << usage;
        return 0;

    }
    if (cmdOptionExists(argv, argv + argc, "-v") || cmdOptionExists(argv, argv + argc, "--version")) {
        // Write version string and exit
        std::cout << my_version << "\n";
	    return 0;
    }

    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::error);
    if (cmdOptionExists(argv, argv + argc, "-d") || cmdOptionExists(argv, argv + argc, "--debug")) {
        logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);
        BOOST_LOG_TRIVIAL(debug) << "Debugging mode engaged.";
    }
    // Default bits per probability is 8, this is the same as the ukbiobank data.
    // 8 bits is probably more accurate from chi anyway
    int number_of_bits = 8;
    if (cmdOptionExists(argv, argv + argc, "-b")) {
        // Pall Melsted knows the maximum error in chi, this needs be roughly as accurate
        // 10 bits do for now, 1/1000 is roughly the chi error (if I understand Arnaldur correctly).
        std::string bitprob = getCmdOption(argv, argv + argc, "-b");
        number_of_bits = std::stoi(bitprob);
    }
    int thread_count = 4;
    if (cmdOptionExists(argv, argv + argc, "-t")) {
        thread_count = std::stoi(getCmdOption(argv, argv + argc, "-t"));
    }
    // Note this is maximum, unreliable, more of a suggestion really.
    omp_set_num_threads(thread_count);
    // Default to 2GB
    // uint64_t mem_use = 2 * 1024 * 1024 * 1024;
    uint64_t mem_use = 2147483648;
    if (cmdOptionExists(argv, argv + argc, "-m")) {
        mem_use = uint64_t(std::stoi(getCmdOption(argv, argv + argc, "-m"))) * 1073741824ll;
    }

    if (cmdOptionExists(argv, argv + argc, "-l") || cmdOptionExists(argv, argv + argc, "--readfile")) {
        std::string infile = getCmdOption(argv, argv + argc, "-l");
        BOOST_LOG_TRIVIAL(debug) << "Reading file " << infile << "\n";
        test_read(infile);
        return 0;
    }
    // Default is to process whole file.
    size_t chunk_size = 0;
    int index = -1;
    if (cmdOptionExists(argv, argv + argc, "-c")) {
        std::string csize = getCmdOption(argv, argv + argc, "-c");
        chunk_size = std::stoi(csize);
        std::string idx = getCmdOption(argv, argv + argc, "-i");
        index = std::stoi(idx);
    }
    std::string infile = getCmdOption(argv, argv + argc, "-f");
    std::string outfile = getCmdOption(argv, argv + argc, "-o");

    // Main thread reader initialized here:
    chiSingleReader creader = chiSingleReader(infile, NULL);
    // Retrieve pn and marker count from chi header.
    size_t marker_count = creader.num_markers();    
    // Retrieve PN count from chi header
    size_t pn_count = creader.num_columns();

    std::vector<std::string> pns;
    read_file(pns, infile + ".pns");

    create_sample_file(outfile, pns);

    std::string indexfname = get_indexfname(infile);

    // Adapt count of markers to process to control parameters.
    // Index refers to chunk number - these are parameters for parallel controller scripts which operate
    // on a higher level, by mass calling of chi2bgen. 
    // So they make little sense if chi2bgen is user-facing
    int main_start = 0;
    
    if (chunk_size != 0) {
        main_start = chunk_size * index;
        if(main_start >= marker_count ) {
            BOOST_LOG_TRIVIAL(debug) << "Index * chunk_size (" << index*chunk_size << ") ge than count of markers.\n";
        }
        marker_count = std::min(chunk_size, marker_count - main_start);
        BOOST_LOG_TRIVIAL(debug) << "Will process markers " << main_start << " to " << main_start + marker_count << "\n";
    }
    int main_end = main_start + marker_count;

    uint64_t chi_marker_size = sizeof(float) * 2 * pn_count;
    uint64_t max_markers_to_buffer = mem_use / chi_marker_size;
    size_t markers_to_buffer = std::min(marker_count, max_markers_to_buffer);

    size_t outer_cycles = (marker_count / markers_to_buffer);
    if (marker_count % markers_to_buffer != 0) {
        outer_cycles += 1;
    }
    
    
    // ------ End of input handling :)
    BOOST_LOG_TRIVIAL(info) << "Converting: " << infile << " to " << outfile << " with version " << my_version << "\n";
    BOOST_LOG_TRIVIAL(info) << "Input chi file has " << pn_count << " pns. Will convert " << marker_count << " markers.\n";


    // Setup some more objects prior to iterating through range.    

    // Prep feedback vars
    size_t report_step = 1000;
    int counter = 0;
    std::chrono::duration<double> seconds;
    auto start = std::chrono::system_clock::now();

    // prep chi input buffers
    size_t chi_inp_size = std::min(marker_count, markers_to_buffer);
    std::vector<float*> chi_inp(chi_inp_size);
    for (size_t i = 0; i < chi_inp.size(); i++) {
        chi_inp[i] = new float[2*pn_count];
    }

    // Indexes must be handled via pointers.
    // Desctructor segfaults also. Just accept memory leak and move on.
    TabixIndex * m_tbx_idx = new TabixIndex(indexfname);
    TabixIndex::Iterator m_tbx_itr = m_tbx_idx->create_iter_range(main_start, main_end);

    // prep bgen variables.
    genfile::bgen::Context mycontext = init_context(pn_count, marker_count);
    uint16_t allele_count = 2;
    genfile::MissingValue missing_allele;
    std::string chr;
    if(std::stoi(m_tbx_itr->field(0)) < 10) {
        chr = "0" + std::string(m_tbx_itr->field(0));
    }
    else {
        chr = std::string(m_tbx_itr->field(0));
    }
    

    // Prep bgen output buffers
    std::vector<std::vector<genfile::byte_t>> mid_buffers(chi_inp_size , std::vector<genfile::byte_t>(0, 0));
    std::vector<std::vector<genfile::byte_t>> gtp_buffers(chi_inp_size, std::vector<genfile::byte_t>(0, 0));
    std::vector<std::vector<genfile::byte_t>> uncompr_gtp_buffers(thread_count, std::vector<genfile::byte_t>(0, 0));

    // prep outfile
    std::ofstream bgenfile;
    bgenfile.open(outfile + ".bgen", std::ios::binary);
    write_bgen_header(bgenfile, mycontext, pns);

    BOOST_LOG_TRIVIAL(debug) << "outer_cycles: " << outer_cycles << "\n";
    BOOST_LOG_TRIVIAL(debug) << "chi_inp_size: " << chi_inp_size << "\n";
    creader.seek(main_start);
    for (size_t k = 0; k < outer_cycles; k++) {

        size_t inner_loop = std::min(marker_count-(k*markers_to_buffer), markers_to_buffer);

        // input loop
        for (size_t i = 0; i < inner_loop; i++) {
            creader.v = chi_inp[i];
            creader.read_next_marker();
        }

        // markerinfo processing loop
        BOOST_LOG_TRIVIAL(debug) << "Iter range of " << main_start + (k*markers_to_buffer) << " to " << main_start + (k*markers_to_buffer) + inner_loop << "\n";
        TabixIndex::Iterator i_tbx_itr = m_tbx_idx->create_iter_range(main_start + (k*markers_to_buffer), main_start + (k*markers_to_buffer) + inner_loop);
        for (size_t i = 0; i < inner_loop; i++) {
            static_cast<void>(genfile::bgen::write_snp_identifying_data(
                &mid_buffers[i],
                mycontext,
                std::string(i_tbx_itr->field(2)),
                std::string(i_tbx_itr->field(2)),
                chr,
                size_t(std::stoi(i_tbx_itr->field(1))),
                2,
                boost::bind(&chi_get_tbx_allele, i_tbx_itr, boost::placeholders::_1)
                ));
                i_tbx_itr++;
        }

        BOOST_LOG_TRIVIAL(debug) << "inner_loop " << inner_loop << "\n";
        // genotype processing loop
        #pragma omp parallel 
        {
            size_t t_inner_loop = inner_loop / omp_get_num_threads();
            size_t t_inner_cursor = omp_get_thread_num() * t_inner_loop;
            if (omp_get_thread_num() == omp_get_num_threads() -1) {
                t_inner_loop += inner_loop % omp_get_num_threads();
            }
            // Temporary container for gtypes until call to writer.finalise()
            // Each thread should have own tbx_itr.
            BOOST_LOG_TRIVIAL(debug) << "t_inner_cursor: " << t_inner_cursor << "\n";
            BOOST_LOG_TRIVIAL(debug) << "t_inner_loop: " << t_inner_loop << "\n";
            for (size_t i = 0; i < t_inner_loop; i++) {            
                float* vals = chi_inp[t_inner_cursor + i];

                genfile::bgen::GenotypeDataBlockWriter mywriter(
                    &uncompr_gtp_buffers[omp_get_thread_num()],
                    &gtp_buffers[t_inner_cursor+i],
                    mycontext,
                    number_of_bits
                );

                mywriter.initialise(pn_count, 2);

                for(size_t j = 0; j < pn_count * 2; j = j + 2) {
                    mywriter.set_sample(j/2);
                    mywriter.set_number_of_entries(2, 4, genfile::ePerPhasedHaplotypePerAllele, genfile::eProbability);

                    if (vals[j] != -1.0) {
                        mywriter.set_value(0, (1 - vals[j]));
                        mywriter.set_value(1, vals[j]);
                    } 
                    else {
                        mywriter.set_value(0, missing_allele);
                        mywriter.set_value(1, missing_allele);
                    }
                    if (vals[j+1] != -1.0) {
                        mywriter.set_value(2, (1 - vals[j+1]));
                        mywriter.set_value(3, vals[j+1]);            
                    }
                    else {
                        mywriter.set_value(2, missing_allele);
                        mywriter.set_value(3, missing_allele);
                    }
                }
                mywriter.finalise();                
                
                // if (++counter % report_step == 0) {
                //     seconds = std::chrono::system_clock::now() - start;
                //     std::cout << counter << " markers were bgen encoded in buffers " << seconds.count() << ".\n";
                // }
            }
        };

        // output loop
        for(size_t i = 0; i < inner_loop; i++) {
            bgenfile.write(reinterpret_cast<char const*>(mid_buffers[i].data()), mid_buffers[i].size() * sizeof(genfile::byte_t));
            bgenfile.write(reinterpret_cast<char const*>(gtp_buffers[i].data()), gtp_buffers[i].size() * sizeof(genfile::byte_t));
        }
    }

    // cleanup
    for ( size_t i = 0; i < chi_inp.size(); i++) {
        delete[] chi_inp[i];
    }

    // delete m_tbx_idx;

    bgenfile.close();
    return 0;
}
