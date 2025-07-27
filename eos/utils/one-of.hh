/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon code which is:
 *     Copyright (c) 2010 Ciaran McCreesh
 *
 * This file is part of the EOS program. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_UTILS_ONE_OF_HH
#define EOS_GUARD_EOS_UTILS_ONE_OF_HH 1

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

namespace eos
{
    /* Helpers */
    struct UnknownTypeForOneOf;

    template <typename Want_, typename... Types_> struct SelectOneOfType;

    template <typename Want_> struct SelectOneOfType<Want_>
    {
            using Type = UnknownTypeForOneOf;
    };

    template <typename Want_, typename Try_, typename... Rest_> struct SelectOneOfType<Want_, Try_, Rest_...>
    {
            using Type = typename std::conditional<std::is_same<Want_, Try_>::value, Try_, typename SelectOneOfType<Want_, Rest_...>::Type>::type;
    };

    /* OneOfVisitor */
    template <typename Type_> struct OneOfVisitorVisit
    {
            virtual void visit(Type_ &) = 0;
    };

    template <typename... Types_> struct OneOfVisitor : OneOfVisitorVisit<Types_>...
    {};

    /* OneOfVisitorWrapper */
    template <typename Visitor_, typename Underlying_, typename Result_, typename... Types_> struct OneOfVisitorWrapperVisit;

    template <typename Visitor_, typename Underlying_, typename Result_> struct OneOfVisitorWrapperVisit<Visitor_, Underlying_, Result_> : Visitor_
    {
            Underlying_ & underlying;

            Result_ result;

            OneOfVisitorWrapperVisit(Underlying_ & u) :
                underlying(u)
            {
            }
    };

    template <typename Visitor_, typename Underlying_> struct OneOfVisitorWrapperVisit<Visitor_, Underlying_, void> : Visitor_
    {
            Underlying_ & underlying;

            OneOfVisitorWrapperVisit(Underlying_ & u) :
                underlying(u)
            {
            }
    };

    template <typename Visitor_, typename Underlying_, typename Type_, typename... Rest_>
    struct OneOfVisitorWrapperVisit<Visitor_, Underlying_, void, Type_, Rest_...> : OneOfVisitorWrapperVisit<Visitor_, Underlying_, void, Rest_...>
    {
            OneOfVisitorWrapperVisit(Underlying_ & u) :
                OneOfVisitorWrapperVisit<Visitor_, Underlying_, void, Rest_...>(u)
            {
            }

            virtual void
            visit(Type_ & t)
            {
                this->underlying.visit(t);
            }
    };

    template <typename Visitor_, typename Underlying_, typename Result_, typename Type_, typename... Rest_>
    struct OneOfVisitorWrapperVisit<Visitor_, Underlying_, Result_, Type_, Rest_...> : OneOfVisitorWrapperVisit<Visitor_, Underlying_, Result_, Rest_...>
    {
            OneOfVisitorWrapperVisit(Underlying_ & u) :
                OneOfVisitorWrapperVisit<Visitor_, Underlying_, Result_, Rest_...>(u)
            {
            }

            virtual void
            visit(Type_ & t)
            {
                this->result = this->underlying.visit(t);
            }
    };

    template <typename Underlying_, typename Result_, typename... Types_>
    struct OneOfVisitorWrapper : OneOfVisitorWrapperVisit<OneOfVisitor<Types_...>, Underlying_, Result_, Types_...>
    {
            OneOfVisitorWrapper(Underlying_ & u) :
                OneOfVisitorWrapperVisit<OneOfVisitor<Types_...>, Underlying_, Result_, Types_...>(u)
            {
            }
    };

    /* OneOf */
    template <typename... Types_> struct OneOfValueBase
    {
            virtual ~OneOfValueBase() = 0;

            virtual void accept(OneOfVisitor<Types_...> &) = 0;
    };

    template <typename... Types_> OneOfValueBase<Types_...>::~OneOfValueBase() = default;

    template <typename Type_, typename... Types_> struct OneOfValue : OneOfValueBase<Types_...>
    {
            Type_ value;

            OneOfValue(const Type_ & type) :
                value(type)
            {
            }

            virtual void
            accept(OneOfVisitor<Types_...> & visitor)
            {
                static_cast<OneOfVisitorVisit<Type_> &>(visitor).visit(value);
            }
    };

    template <typename... Types_> class OneOf
    {
        private:
            std::shared_ptr<OneOfValueBase<Types_...>> _value;

            OneOf(OneOfValueBase<Types_...> * ptr) :
                _value(ptr)
            {
            }

        public:
            OneOf() {}

            template <typename Type_>
            OneOf(const Type_ & value) :
                _value(new OneOfValue<typename SelectOneOfType<Type_, Types_...>::Type, Types_...>{ value })
            {
            }

            OneOf(const OneOf & other) :
                _value(other._value)
            {
            }

            bool
            empty() const
            {
                return ! _value;
            }

            template <typename Type_>
            OneOf &
            operator= (const Type_ & value)
            {
                _value.reset(new OneOfValue<typename SelectOneOfType<Type_, Types_...>::Type, Types_...>{ value });
                return *this;
            }

            OneOf &
            operator= (const OneOf & other)
            {
                _value = other._value;

                return *this;
            }

            template <typename Visitor_>
            void
            accept(Visitor_ & visitor) const
            {
                OneOfVisitorWrapper<Visitor_, void, Types_...> visitor_wrapper(visitor);

                _value->accept(visitor_wrapper);
            }

            template <typename Result_, typename Visitor_>
            Result_
            accept_returning(Visitor_ & visitor) const
            {
                OneOfVisitorWrapper<Visitor_, Result_, Types_...> visitor_wrapper(visitor);

                _value->accept(visitor_wrapper);

                return visitor_wrapper.result;
            }
    };
} // namespace eos

#endif
