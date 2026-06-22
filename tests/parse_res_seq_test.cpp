/**
 * @file parse_res_seq_test.cpp
 * @brief Unit tests for scoring::parse_res_seq.
 */
#include <string>

#include <pairfinder/scoring/torsions.hpp>

#include "test_helpers.hpp"

namespace {

using pairfinder::scoring::parse_res_seq;
using pairfinder::tests::check;
using pairfinder::tests::finish;

}  // namespace

int main() {
    {
        const auto p = parse_res_seq("A-G-20");
        check(p.has_value(), "A-G-20: parsed");
        if (p) {
            check(p->chain == "A", "A-G-20: chain");
            check(p->num == 20, "A-G-20: num");
            check(p->icode.empty(), "A-G-20: icode");
        }
    }

    {
        const auto p = parse_res_seq("G-C--8");
        check(p.has_value(), "G-C--8: parsed");
        if (p) {
            check(p->chain == "G", "G-C--8: chain");
            check(p->num == -8, "G-C--8: num");
        }
    }

    {
        const auto p = parse_res_seq("A-G-20A");
        check(p.has_value(), "A-G-20A: parsed");
        if (p) {
            check(p->chain == "A", "A-G-20A: chain");
            check(p->num == 20, "A-G-20A: num");
            check(p->icode == "A", "A-G-20A: icode");
        }
    }

    {
        const auto p = parse_res_seq("B-DA-5");
        check(p.has_value(), "B-DA-5: parsed");
        if (p) {
            check(p->chain == "B", "B-DA-5: chain");
            check(p->num == 5, "B-DA-5: num");
        }
    }

    {
        check(!parse_res_seq("malformed").has_value(), "malformed: nullopt");
        check(!parse_res_seq("A-G-").has_value(), "missing seq: nullopt");
        check(!parse_res_seq("no-dash").has_value(), "no dash: nullopt");
    }

    return finish("parse_res_seq_test");
}
