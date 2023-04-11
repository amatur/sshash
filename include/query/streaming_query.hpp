#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"

#include "../gz/zip_stream.hpp"
#include "streaming_query_canonical_parsing.hpp"
#include "streaming_query_regular_parsing.hpp"

namespace sshash {

template <typename Query>
void streaming_query_from_fasta_file_multiline(dictionary const* dict, std::istream& is) {
    buffered_lines_iterator it(is);
    std::string buffer;
    uint64_t k = dict->k();
    Query query(dict);
    query.start();
    while (!it.eof()) {
        bool empty_line_was_read = it.fill_buffer(buffer);
        for (uint64_t i = 0; i != buffer.size() - k + 1; ++i) {
            char const* kmer = buffer.data() + i;
            auto answer = query.lookup_advanced(kmer);
            if( answer.kmer_id != constants::invalid_uint64){
                std::cout<<answer.kmer_id<<std::endl;
            }else{
                std::cout<<-1<<std::endl;
            }
        }
        if (empty_line_was_read) { /* re-start the kmers' buffer */
            buffer.clear();
            query.start();
        } else {
            if (buffer.size() > k - 1) {
                std::copy(buffer.data() + buffer.size() - k + 1, buffer.data() + buffer.size(),
                          buffer.data());
                buffer.resize(k - 1);
            }
        }
    }
    return;
}


void dictionary::streaming_query_from_file(std::string const& filename) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_from_fasta_file_multiline<streaming_query_canonical_parsing>(this, is);
    is.close();
    return;
}

}  // namespace sshash