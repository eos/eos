/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/statistics/chain-group.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/verify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace eos
{
    template <> struct Implementation<ChainGroup>
    {
        ChainGroup::RValueFunction rvalue_function;

        double max_rvalue;

        unsigned number_of_parameters;

        std::vector<HistoryPtr> chains;

        // indices of chains which are added
        std::vector<unsigned> chain_indices;

        // indices of parameters whose R-value ought to be checked
        std::vector<unsigned> parameter_indices;

        // one vector of parameter means/variances for each chain
        std::vector<std::vector<double>> parameter_means;
        std::vector<std::vector<double>> parameter_variances;

        VerifiedRange<double> skip_initial;

        Implementation(const ChainGroup::RValueFunction & rvalue_function, const double & max_rvalue,
                       const HistoryPtr & initial_chain, const unsigned & index, const double & skip_initial) :
            rvalue_function(rvalue_function),
            max_rvalue(max_rvalue),
            number_of_parameters(initial_chain->states.front().point.size()),
            parameter_indices(number_of_parameters),
            skip_initial(0, 1, skip_initial)
        {
            std::iota(parameter_indices.begin(), parameter_indices.end(), 0);
            add(initial_chain, index);
        }

        ~Implementation()
        {
        }

        void add(const HistoryPtr & chain, const unsigned & index)
        {
            chains.push_back(chain);
            chain_indices.push_back(index);

            // calculate mean
            std::vector<double> means, variances;
            MarkovChain::State::Iterator it_begin = chain->states.cbegin() + unsigned(skip_initial * chain->states.size());
            MarkovChain::State::Iterator it_end = chain->states.cend();
            chain->mean_and_variance(it_begin, it_end, means, variances);

            parameter_means.push_back(means);
            parameter_variances.push_back(variances);
        }

        bool overlaps(const HistoryPtr & chain) const
        {
            std::vector<double> all_chain_means(chains.size() + 1, 0.0), all_chain_variances(chains.size() + 1, 0.0);

            // todo should be another type of Exception
            if (chain->states.front().point.size() != number_of_parameters)
                throw InternalError("cluster: chain size doesn't match");

            // compute statistics for the chain to test
            std::vector<double> new_chain_means(number_of_parameters, 0.0);
            std::vector<double> new_chain_variances(new_chain_means);
            MarkovChain::State::Iterator it_begin = chain->states.cbegin() + unsigned(skip_initial * chain->states.size());
            chain->mean_and_variance(it_begin, chain->states.cend(), new_chain_means, new_chain_variances);

            // suppose n=10 and skip = 15%, then in the iterators, only one elements is skipped,
            // but here the length would be 8 instead of 9, so use ceiling
            unsigned number_of_points = std::ceil((1.0 - skip_initial) * chain->states.size());

            // check overlap in each parameter dimension
            for (auto i = parameter_indices.cbegin() ; i != parameter_indices.cend() ; ++i)
            {
                all_chain_means.clear();
                all_chain_variances.clear();

                // consider the all_chain_means/all_chain_variances of all chains already in the cluster
                unsigned chain_index = 0;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++chain_index)
                {
                    all_chain_means.push_back(parameter_means[chain_index][*i]);
                    all_chain_variances.push_back(parameter_variances[chain_index][*i]);
                }

                // and compare with the new chain
                all_chain_means.push_back(new_chain_means[*i]);
                all_chain_variances.push_back(new_chain_variances[*i]);

                double rvalue = rvalue_function(all_chain_means, all_chain_variances, number_of_points);

                if (rvalue > max_rvalue)
                {
                    Log::instance()->message("Cluster.overlaps", ll_debug)
                        << "Parameter " << *i << ": r value too large (" << rvalue << " > " << max_rvalue << ")";
                    return false;
                }
            }

            return true;
        }
    };

    ChainGroup::ChainGroup(const RValueFunction & rvalue_function, const double & max_rvalue, const HistoryPtr & initial_chain, const unsigned & index, const double & skip_initial) :
        PrivateImplementationPattern<ChainGroup>(new Implementation<ChainGroup>(rvalue_function, max_rvalue, initial_chain, index, skip_initial))
    {
    }

    ChainGroup::~ChainGroup()
    {
    }

    void
    ChainGroup::add(const HistoryPtr & chain, const unsigned & index)
    {
        _imp->add(chain, index);
    }

    ChainGroup::Iterator
    ChainGroup::begin() const
    {
        return ChainGroup::Iterator(_imp->chains.begin());
    }

    ChainGroup::Iterator
    ChainGroup::end() const
    {
        return ChainGroup::Iterator(_imp->chains.end());
    }

    ChainGroup::IndexIterator
    ChainGroup::begin_indices() const
    {
        return ChainGroup::IndexIterator(_imp->chain_indices.begin());
    }

    ChainGroup::IndexIterator
    ChainGroup::end_indices() const
    {
        return ChainGroup::IndexIterator(_imp->chain_indices.end());
    }

    bool
    ChainGroup::overlaps(const HistoryPtr & chain) const
    {
        return _imp->overlaps(chain);
    }

    void
    ChainGroup::parameter_indices(const std::vector<unsigned> & indices)
    {
        // remove old indices
        _imp->parameter_indices.clear();

        // add new indices in ascending order
        std::vector<unsigned> indices_sorted(indices.cbegin(), indices.cend());
        std::sort(indices_sorted.begin(), indices_sorted.end());

        for(auto i = indices_sorted.cbegin() ; i != indices_sorted.cend() ; ++i)
        {
            if (*i > _imp->number_of_parameters)
                throw InternalError("Cluster::parameter_indices: index " + stringify(*i)
                                    + " out of range");
            _imp->parameter_indices.push_back(*i);
        }
    }

    std::vector<double>
    ChainGroup::mean() const
    {
        std::vector<double> result(_imp->parameter_means.front());

        // build average using Welford's method
        double n = 1;
        for (auto par_mean = _imp->parameter_means.cbegin() + 1 ; par_mean != _imp->parameter_means.cend() ; ++par_mean, ++n)
        {
            for (unsigned i = 0 ; i < par_mean->size() ; ++i)
            {
                result[i] += ((*par_mean)[i] - result[i]) / n;
            }
        }
        return result;
    }

    const std::vector<std::vector<double>> &
    ChainGroup::means() const
    {
        return _imp->parameter_means;
    }

    const std::vector<std::vector<double>> &
    ChainGroup::variances() const
    {
        return _imp->parameter_variances;
    }

    template <>
    struct WrappedForwardIteratorTraits<ChainGroup::IteratorTag>
    {
        typedef std::vector<HistoryPtr>::iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ChainGroup::IteratorTag, HistoryPtr>;

    template <>
    struct WrappedForwardIteratorTraits<ChainGroup::IndexIteratorTag>
    {
        typedef std::vector<unsigned>::iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ChainGroup::IndexIteratorTag, unsigned>;
}
