#pragma once
#ifndef BMSRS_SORT_H
#define BMSRS_SORT_H
#include <iostream>
#include <cassert>
#include <bit>
#include <cstdint>

namespace Algo
{ enum : uint32_t
{
    InsertionThreshold = 88,
    AscendingThreshold = 8,
    LargeDataThreshold = 128,
    BlockSize          = 64, 
#if __cpp_lib_bitops >= 201907L
    DoubleWordBitCount = 31,
#else
    DeBruijnShiftAmoun = 58
#endif
};

/**
 * The DeBruijn constant.
 */
#if __cpp_lib_bitops < 201907L
constexpr uint64_t DeBruijn64 =
    0x03F79D71B4CB0A89L;
#endif

/**
 * The DeBruijn map from key to integer
 * square index.
 */
#if __cpp_lib_bitops < 201907L
constexpr uint8_t DeBruijnTableF[] = 
{
    0,  47,  1, 56, 48, 27,  2, 60,
    57, 49, 41, 37, 28, 16,  3, 61,
    54, 58, 35, 52, 50, 42, 21, 44,
    38, 32, 29, 23, 17, 11,  4, 62,
    46, 55, 26, 59, 40, 36, 15, 53,
    34, 51, 20, 43, 31, 22, 10, 45,
    25, 39, 14, 33, 19, 30,  9, 24,
    13, 18,  8, 12,  7,  6,  5, 63
};
#endif

/**
 * Fill trailing bits using prefix fill.
 *
 * @code
 * Example:
 *       10000000 >> 1
 *     = 01000000 | 10000000
 *     = 11000000 >> 2
 *     = 00110000 | 11000000
 *     = 11110000 >> 4
 *     = 00001111 | 11110000
 *     = 11111111
 * @endcode
 * @tparam E The type
 * @param x The integer
 */
#if __cpp_lib_bitops < 201907L
constexpr void parallelPrefixFill
    (
    uint32_t & x
    ) 
{
    x |= x >> 1U;
    x |= x >> 2U;
    x |= x >> 4U;
    x |= x >> 8U;
    x |= x >> 16U;
}
#endif

/**
 * Calculates floor of log2
 * 
 * @authors Kim Walisch - source
 * @authors Mark Dickinson - source
 * @authors Ellie Moore 
 * @param bb bitboard to scan
 * @precondition bb != 0
 * @return index (0..63) of most significant one bit
 */
constexpr int log2
    (
    uint32_t l
    ) 
{
    assert(l != 0);
#if __cpp_lib_bitops >= 201907L
    return std::countl_zero(l) ^ DoubleWordBitCount; 
#else
    parallelPrefixFill(l);
    return DeBruijnTableF[(int)
        ((l * DeBruijn64) >> DeBruijnShiftAmoun)
    ];
#endif
}

struct interval final { void * a, * b; };

/**
 * A simple swap method.
 *
 * @tparam E the element type
 * @param i the first element pointer
 * @param j the second element pointer
 */
template<typename E>
constexpr void swap
    (
    E *const i,
    E *const j
    ) 
{
    E const el = *i; *i = *j; *j = el;
}

/**
 * A generic "sift down" method (AKA max-heapify)
 *
 * @tparam E the element type
 * @param a the pointer to the base of the current
 * sub-array
 * @param i the starting index
 * @param size the size of the current sub-array
 */
template<typename E, class Cmp>
inline void siftDown
    (
    E* const a,
    const int i,
    const int size,
    const Cmp cmp
    ) 
{
    // Store size in
    // a local variable.
    const size_t n = size;

    // Establish non-leaf
    // boundary.
    const size_t o = n >> 1U;

    // Extract the element
    // to sift.
    E z = a[i];

    // initialize temporary
    // variables.
    size_t x = i, l, r;

    // consider only non-leaf
    // nodes.
    while(x < o) 
    {
        // y is currently
        // left child element.
        // Note: "l" here stands
        // for "left"
        r = (l = (x << 1U) + 1) + 1;
        E y = a[l];

        // if right child is
        // within the heap...
        // AND
        // if right child element
        // is greater than left
        // child element,
        // THEN
        // assign right child to
        // y and right index to l.
        // Note: "l" now stands
        // for "larger"
        if(r < n && cmp(y, a[r]))
            y = a[l = r];

        // if y is less than or
        // equal to the element
        // we are sifting, then
        // we are done.
        if(!cmp(z, y)) break;

        // move y up to the
        // parent index.
        a[x] = y;

        // Set parent index to
        // be the index of
        // the largest child.
        x = l;
    }

    // Place the sifted element.
    a[x] = z;
}

/**
 * <h1>
 *  <b>
 *  <i>Heap Sort</i>
 *  </b>
 * </h1>
 *
 * <p>
 * Classical heap sort that sorts the given range
 * in ascending order, building a max heap and
 * continuously sifting/swapping the max element
 * to the previous rightmost index.
 * </p>
 * @author Ellie Moore
 * @tparam E the element type
 * @param low a pointer to the leftmost index
 * @param high a pointer to the rightmost index
 */
template<typename E, class Cmp>
inline void hSort
    (
    E* const low,
    E* const high,
    const Cmp cmp
    ) 
{
    E* r = high + 1;
    E* const l = low;

    // Build the heap.
    int x = r - l;
    for(int i =
        (x >> 1U); i >= 0; --i)
        siftDown(l, i, x, cmp);
    
    // Sort.
    while(l < --r) 
    {
        const E z = *l; *l = *r;
        siftDown(l, 0, --x, cmp);
        *r = z;
    }
}

/**
 * <h1>
 *  <b>
 *  <i>Insertion Sort</i>
 *  </b>
 * </h1>
 *
 * <p>
 * Classical ascending insertion sort packaged with a
 * "pairing" optimization to be used in the context of
 * Quicksort.
 * </p>
 * 
 * <p>
 * This optimization is used whenever the portion of
 * the array to be sorted is padded on the left by
 * a portion with lesser elements. The fact that all of
 * the elements on the left are automatically less than
 * the elements in the current portion allows us to skip
 * the costly lower boundary check in the nested loops
 * and insert two elements in one go.
 * </p>
 *
 * @authors Josh Bloch - source
 * @authors Jon Bently - source
 * @authors Orson Peters - source
 * @authors Ellie Moore
 * @tparam E the element type
 * @tparam Are we sorting optimistically?
 * @param leftmost whether this is the leftmost part
 * @param low a pointer to the leftmost index
 * @param high a pointer to the rightmost index
 * left-most partition.
 */
template
<bool NoGuard, bool Guard, bool Bail = true, typename E, class Cmp>
inline bool iSort
    (
    E *const low, 
    E *const high,
    const Cmp cmp,
    const bool leftmost = true
    ) 
{
    E* l = low;
    E* r = high;
    int moves = 0;

    // We aren't guarding, jump
    // straight into pair insertion
    // sort.
    if constexpr (Guard)
        goto g1;
    if constexpr (NoGuard)
        goto g2;

    if (leftmost) 
    {
        g1:

        // Traditional
        // insertion
        // sort.
        for (E *i = l + 1; i <= r; ++i) 
        {
            E t = *i, *j = i - 1;
            for (; j >= l && cmp(t, *j); --j)
                j[1] = *j;
            j[1] = t;

            if constexpr (Bail) 
            {
                // If we have moved too
                // many elements, abort.
                moves += (i - 1) - j;
                if(moves > AscendingThreshold)
                    return false;
            }
        }
    } 
    else 
    {
        g2:

        // Pair insertion sort.
        // Skip elements that are
        // in ascending order.
        do if (l++ >= r) return true;
        while (!cmp(*l, *(l - 1)));

        // This sort uses the sub
        // array at left to avoid
        // the lower bound check.
        // Assumes that this is not
        // the leftmost partition.
        for (E *i = l; ++l <= r; i = ++l) 
        {
            E ex = *i, ey = *l;

            // Make sure that
            // we insert the
            // larger element
            // first.
            if (cmp(ey, ex)) 
            {
                ex = ey;
                ey = *i;
                if constexpr (Bail) 
                    ++moves;
            }

            // Insert the two
            // in one downward
            // motion.
            while (cmp(ey, *--i))
                i[2] = *i;
            (++i)[1] = ey;
            while (cmp(ex, *--i))
                i[1] = *i;
            i[1] = ex;

            if constexpr (Bail) 
            {
                // If we have moved too
                // many elements, abort.
                moves += (l - 2) - i;
                if(moves > AscendingThreshold) 
                    return false;
            }
        }

        // For odd length arrays,
        // insert the last element.
        E ez = *r;
        while (cmp(ez, *--r))
            r[1] = *r;
        r[1] = ez;
    }
    return true;
}

/**
 * Explicit constexpr ternary.
 * 
 * @tparam EXP the constexpr
 * @tparam E the return type
 * @param a the true value
 * @param b the false value
 */
template<bool EXP, typename E>
constexpr E ternary
    (
    E a, 
    E b
    )
{
    if constexpr (EXP) return a;
    else return b;
}

/**
 * Scramble a few elements to help
 * break patterns.
 *
 * @tparam E the element type
 * @param i the first element pointer
 * @param j the second element pointer
 */
template<typename E>
constexpr void scramble
    (
    E* const low, 
    E* const high, 
    const size_t len
    ) 
{
    if(len >= InsertionThreshold) 
    {
        const int _4th = len >> 2U;
        swap(low,  low  + _4th);
        swap(high, high - _4th);
        if(len > LargeDataThreshold)
        {
            swap(low  + 1, low  + (_4th + 1));
            swap(low  + 2, low  + (_4th + 2));
            swap(high - 2, high - (_4th + 2));
            swap(high - 1, high - (_4th + 1));
        }
    }
}

/**
 * Aligns the given pointer on 64-bit 
 * boundary.
 * 
 * @tparam E the element type
 * @param p pointer to memory to align
 */
template<typename E>
constexpr E* align
    (
    E* p
    )
{
    return reinterpret_cast<E*>((
        reinterpret_cast<uintptr_t>(p) + (BlockSize - 1)
    ) & -uintptr_t(BlockSize));
}

/**
 * The sort stack(s). 
 */
interval stack[log2(UINT32_MAX)];
uint8_t heights[log2(UINT32_MAX)];
interval * end = stack + log2(UINT32_MAX);

/**
 * <h1>
 *  <b>
 *  <i>Blipsort</i>
 *  </b>
 * </h1>
 *
 * <h2>Branchless Lomuto</h2>
 * <p>
 * The decades-old partitioning algorithm recently 
 * made a resurgence when researchers discovered 
 * ways to remove the inner branch. Lukas Bergdoll
 * and Orson Peters' method— published a little 
 * under two months ago— is the fastest yet. It 
 * employs a gap in the data to move elements 
 * twice per iteration rather than swapping them 
 * (three moves). For arithmetic and pointer types, 
 * Blipsort employs branchless Lomuto partitioning. 
 * For other, larger types, Blipsort uses branchless 
 * or branchful Hoare partitioning.
 * </p>
 * 
 * <h2>Pivot Selectivity</h2>
 * <p>
 * Blipsort carefully selects the pivot from the 
 * middle of five sorted candidates. These 
 * candidates allow the sort to determine whether 
 * the data in the current interval is approximately 
 * descending and inform its "partition left" strategy.
 * </p>
 * 
 * <h2>Insertion Sort</h2>
 * <p>
 * Blipsort uses Insertion sort on small intervals 
 * where asymptotic complexity matters less and 
 * instruction overhead matters more. Blipsort 
 * employs Java's Pair Insertion sort on every 
 * interval except the leftmost. Pair insertion 
 * sort inserts two elements at a time and doesn't 
 * need to perform a lower bound check, making it 
 * slightly faster than normal insertion sort in 
 * the context of quicksort.
 * </p>
 * 
 * <h2>Pivot Retention</h2>
 * <p>
 * Similar to PDQsort, if any of the three middlemost 
 * candidate pivots is equal to the rightmost element 
 * of the partition at left, Blipsort moves equal 
 * elements to the left with branchless Lomuto and 
 * continues to the right, solving the dutch-flag 
 * problem and yeilding linear time on data comprised 
 * of equal elements.
 * </p> 
 * 
 * <h2>Optimism</h2>
 * <p>
 * Similar to PDQsort, if the partition is "good" 
 * (not highly unbalanced), Blipsort switches to 
 * insertion sort. If the Insertion sort makes more 
 * than a constant number of moves, Blipsort bails 
 * and resumes quicksort. This allows Blipsort to 
 * achieve linear time on already-sorted data.
 * </p>
 * 
 * <h2>Breaking Patterns</h2>
 * <p>
 * Like PDQsort, if the partition is bad, Blipsort 
 * scrambles some elements to break up patterns.
 * </p>
 * 
 * <h2>Rotation</h2>
 * <p>
 * When all of the candidate pivots are strictly 
 * descending, it is very likely that the interval 
 * is descending as well. Lomuto partitioning slows 
 * significantly on descending data. Therefore, 
 * Blipsort neglects to sort descending candidates 
 * and instead swap-rotates the entire interval 
 * before partitioning.
 * </p>
 *
 * <h2>Custom Comparators</h2>
 * <p>
 * Blipsort allows its user to implement a custom 
 * boolean comparator. A comparator is best 
 * implemented with a lambda and no branches. A 
 * comparator implemented with a lambda can be 
 * inlined by an optimizing compiler, while a 
 * comparator implemented with a constexpr/inline 
 * function typically cannot.
 * </p>
 *
 * @authors Josh Bloch - source
 * @authors Jon Bently - source
 * @authors Orson Peters - source
 * @authors Lukas Bergdoll - source
 * @authors Stefan Edelkamp - source
 * @authors Armin Weiß - source
 * @authors Ellie Moore
 * @tparam E the element type
 * @tparam Root whether this is the sort root
 * @param leftmost whether this is the leftmost part
 * @param low a pointer to the leftmost index
 * @param high a pointer to the rightmost index
 * @param height the distance of the current sort
 * tree from the initial height of 2log<sub>2</sub>n
 */
template
<bool Expense, bool Block, bool Root = true, typename E, class Cmp>
void qSort
    (
    E * low,
    E * high,
    int height,
    const Cmp cmp
    ) 
{
    // Set up the stack(s).
    interval * s = stack;
    uint8_t * sh = heights;

    // Save low.
    E* const olow = low;

    // Tail call loop.
    for(uint32_t x = high - low;;) 
    {
        // Declare variables
        // ahead of time to
        // avoid issues with goto.
        E * l, * k, * g, p;
        uint32_t ls, gs, _8th;
        bool work;

        // Find an inexpensive
        // approximation of a third of
        // the interval.
        const uint32_t y = x >> 2U,
            _3rd = y + (y >> 1U),
            _6th = _3rd >> 1U;

        // Find an approximate
        // midpoint of the interval.
        E *const mid = low + (x >> 1U);

        // Assign tercile indices
        // to candidate pivots.
        E *const sl = low  + _3rd;
        E *const sr = high - _3rd;

        // Assign outer indices
        // to candidate pivots.
        E * cl = low  + _6th;
        E * cr = high - _6th;

       // If the candidates aren't
        // descending...
        // Insertion sort all five
        // candidate pivots in-place.
        if((!cmp(*cl, *low)) | 
           (!cmp(*sl, *cl))  |
           (!cmp(*mid, *sl)) |
           (!cmp(*sr, *mid)) |
           (!cmp(*cr, *sr))  |
           (!cmp(*high, *cr)))
        {
            
            if(cmp(*low, *cl)) 
                cl = low;
            if(cmp(*cr, *high)) 
                cr = high;

            if (cmp(*sl, *cl)) 
            {
                E e = *sl;
                *sl = *cl;
                *cl =   e;
            }

            if (cmp(*mid, *sl)) 
            {
                E e  = *mid;
                *mid =  *sl;
                *sl  =    e;
                if (cmp(e, *cl)) 
                {
                    *sl = *cl;
                    *cl =   e;
                }
            }

            if (cmp(*sr, *mid))
            {
                E e  =  *sr;
                *sr  = *mid;
                *mid =    e;
                if (cmp(e, *sl))
                {
                    *mid = *sl;
                    *sl  =   e;
                    if (cmp(e, *cl)) 
                    {
                        *sl = *cl;
                        *cl =   e;
                    }
                }
            }

            if (cmp(*cr, *sr))
            {
                E e = *cr;
                *cr = *sr;
                *sr =   e;
                if (cmp(e, *mid)) 
                {
                    *sr  = *mid;
                    *mid =    e;
                    if (cmp(e, *sl)) 
                    {
                        *mid = *sl;
                        *sl  =   e;
                        if (cmp(e, *cl)) 
                        {
                            *sl = *cl;
                            *cl =   e;
                        }
                    }
                }
            }
        }

        // If the candidates are
        // descending, then the
        // interval is likely to
        // be descending somewhat.
        // rotate the entire interval
        // around the midpoint.
        // Don't worry about the
        // even size case. One
        // out-of-order element
        // is no big deal.
        else
        {
            E* u = low;
            E* q = high;
            while(u < mid)
            {
                E e = *u;
                *u++ = *q;
                *q-- = e;
            }
        }
        
        // If any middle candidate 
        // pivot is equal to the 
        // rightmost element of the 
        // partition to the left,
        // swap pivot duplicates to 
        // the side and sort the 
        // remainder. This is an
        // alternative to dutch-flag
        // partitioning.
        if(low != olow)
        {
            // Check the pivot to 
            // the left.
            E h = *(low - 1);
            if(!cmp(h, *sl)  || 
               !cmp(h, *mid) || 
               !cmp(h, *sr)) 
            {
               E* l = low - 1,
                 * g = high + 1;

                // skip over data
                // in place.         
                while(cmp(h, *--g));

                if(g == high)
                    while(!cmp(h, *++l) && l < g);
                else 
                    while(!cmp(h, *++l));

                // If we are sorting 
                // non-arithmetic types,
                // use Hoare for fewer
                // moves.
                if constexpr (Expense)
                {
        /**
         * Partition left by branchful Hoare scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |    ... == h     |        ... ? ...        |     ... > h     |
         * +-------------------------------------------------------------+
         * ^                 ^                         ^                 ^
         * low               l                         k              high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... == h           |            > h ...           |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
                    while(l < g)
                    {
                        swap(l, g);
                        while(cmp(h, *--g));
                        while(!cmp(h, *++l));
                    }
                }

                // If we are sorting 
                // arithmetic types,
                // use branchless lomuto
                // for fewer branches.
                else
                {   
        /**
         * Partition left by branchless Lomuto scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |  ... == h  |  ... > h  | * |     ... ? ...      |  ... > h  |
         * +-------------------------------------------------------------+
         * ^            ^           ^                        ^           ^
         * low          l           k                        g         high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... == h           |            > h ...           |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
                    E * k = l, p = *l;
                    while(k < g)
                    {
                        *k = *l;
                        *l = *++k;
                        l += !cmp(h, *l);
                    }
                    *k = *l;
                    *l = p;
                    l += !cmp(h, p);
                }

                // Advance low to the
                // start of the right
                // partition.               
                low = l;

                // If we have nothing 
                // left to sort, pop
                // the next interval
                // from the stack.
                if(low >= high)
                    goto l3;

                // Calculate the interval 
                // width and loop.
                x = high - low;

                // If the interval is small,
                // insertion sort.
                if(x < InsertionThreshold)
                {
                    iSort<1,0,0>
                    (low, high, cmp);
                    goto l3;
                }

                // Continue with quicksort.
                continue;
            }
        }

        // Initialize l and k.
        l = ternary<Expense>(low, low - 1);
        k = high + 1;

        // Assign midpoint to pivot
        // variable.
        p = *mid;

        // If we are sorting 
        // non-arithmetic types, bring 
        // left end inside. Left end 
        // will be replaced and pivot 
        // will be swapped back later.
        if constexpr(Expense)
            *mid = *low;

        // skip over data
        // in place.
        while(cmp(*++l, p));

        // If we are sorting 
        // arithmetic types, bring 
        // left end inside. Left end 
        // will be replaced and pivot 
        // will be swapped back later.
        if constexpr(!Expense)
            *mid = *l;

        // skip over data
        // in place.
        if(ternary<Expense>
            (l == low + 1, l == low))
            while(!cmp(*--k, p) && k > l);
        else 
            while(!cmp(*--k, p));

        // Will we do a significant 
        // amount of work during 
        // partitioning?
        work = 
        ((l - low) + (high - k)) 
            < (x >> 1U);

        // If we are sorting 
        // non-arithmetic types and
        // conserving memory, use 
        // Hoare for fewer moves.
        if constexpr (Expense && !Block)
        {
        /**
         * Partition by branchful Hoare scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |     ... < p     |        ... ? ...        |    ... >= p     |
         * +-------------------------------------------------------------+
         * ^                 ^                         ^                 ^
         * low               l                         k              high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... < p            |            >= p ...          |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
            while(l < k)
            {
                swap(l, k);
                while(cmp(*++l, p));
                while(!cmp(*--k, p));
            }
            *low = *--l; *l = p;
        }

        // If we are sorting 
        // non-arithmetic types and
        // not conserving memory, use 
        // Block Hoare for fewer moves
        // and fewer branches.
        else if constexpr (Expense)
        {
        /**
         * Partition by branchless (Block) Hoare scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |  ... < p  |   cmp   |   ... ? ...   |   cmp   |  ... >= p   |
         * +-------------------------------------------------------------+
         * ^           ^         ^               ^         ^             ^
         * low         _low      l               k     _high          high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... < p            |            >= p ...          |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
            if(l < k)
            {
                swap(l++, k);

                // Set up blocks and
                // align base pointers to
                // the cacheline.
                uint8_t 
                ols[BlockSize << 1U], 
                oks[BlockSize << 1U]; 
                uint8_t
                * olp = align(ols), 
                * okp = align(oks); 
                
                // Initialize frame pointers.
                E * _low = l, * _high = k;

                // Initialize offset counts and
                // start indices for swap routine.
                size_t nl = 0, nk = 0, ls = 0, ks = 0;
                
                while(l < k) 
                {
                    // If both blocks are empty, split 
                    // the interval in two. Otherwise
                    // give the whole interval to one
                    // block.
                    size_t xx = k - l,
                    lspl = -(nl == 0) & (xx >> (nk == 0)),
                    kspl = -(nk == 0) & (xx - lspl);
                    
                    // Fill the offset blocks. If the split 
                    // for either block is larger than 64,
                    // crop it and unroll the loop. Otherwise,
                    // keep the loop fully rolled. This should
                    // only happen near the end of partitioning.
                    if(lspl >= BlockSize)
                    {
                        size_t i = -1;
                        do
                        {
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                            olp[nl] = ++i; nl += !cmp(*l++, p);
                        } while(i < BlockSize - 1);
                    }
                    else
                        for(size_t i = 0; i < lspl; ++i)
                            olp[nl] = i, nl += !cmp(*l++, p);

                    if(kspl >= BlockSize)
                    {
                        size_t i = 0;
                        do
                        {
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                            okp[nk] = ++i; nk += cmp(*--k, p);
                        } while(i < BlockSize);

                    }
                    else
                        for(size_t i = 0; i < kspl;)
                            okp[nk] = ++i, nk += cmp(*--k, p);

                    // n = min(nl, nk), branchless.
                    size_t n = 
                        (nl & -(nl < nk)) + (nk & -(nl >= nk));

                    // Swap the elements using the offsets. 
                    // Set up working block pointers and lower 
                    // block end pointer.
                    uint8_t* 
                    ll = olp + ls, * kk = okp + ks, * e = ll + n;

                    // If the offset counts are equal, we are likely
                    // to be ascending or descending. If ascending,
                    // we don't need to do anything. If descending, 
                    // use swaps to stay O(n). Both blocks must 
                    // contain n offsets. If either block is empty,
                    // fill it and come back.
                    if(nl == nk)
                        for(; ll < e; ++ll, ++kk)
                            swap(_low + *ll, _high - *kk);

                    // Otherwise, swap using a cyclic permutation.
                    // Both blocks must contain n offsets. If either 
                    // block is empty, fill it and come back.
                    else if(n > 0)
                    {
                        E* _l = _low + *ll, * _k = _high - *kk;
                        E t = *_l; *_l = *_k;
                        for(++ll, ++kk; ll < e; ++ll, ++kk)
                        {
                            _l = _low  + *ll; *_k = *_l;
                            _k = _high - *kk; *_l = *_k;
                        } 
                        *_k = t;
                    }

                    // Adjust offset counts and starts. If a block
                    // is empty, adjust its frame pointer.
                    nl -= n; nk -= n;
                    if(nl == 0) { ls = 0; _low  = l; } else ls += n; 
                    if(nk == 0) { ks = 0; _high = k; } else ks += n;
                }

                // swap the remaining elements into place.
                if(nl)
                {
                    olp += ls;
                    for(uint8_t* ll = olp + nl;;)
                    {
                        swap(_low + *--ll, --l); 
                        if(ll <= olp) break;
                    }
                }

                if(nk)
                {
                    okp += ks;
                    for(uint8_t* kk = okp + nk;;)
                    {
                        swap(_high - *--kk, l++);
                        if(kk <= okp) break;
                    }
                }
            }
            
            // Move the pivot into place.
            *low = *--l; *l = p;
        }

        // If we are sorting 
        // arithmetic types, use 
        // branchless lomuto for 
        // fewer branches.
        else
        {
        /**
         * Partition by branchless Lomuto scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |  ... < p  |  ... >= p  | * |     ... ? ...     |  ... >= p  |
         * +-------------------------------------------------------------+
         * ^           ^            ^                       ^            ^
         * low         l            g                       k         high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... < p            |            >= p ...          |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
            g = l; 
            
            // If we are not conserving 
            // memory, unroll the
            // loop for a tiny boost.
            if constexpr (Block)
            {
                E* u = k - (BlockSize >> 2U);
                while(g < u)
                {
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                    *g = *l; *l = *++g; l += cmp(*l, p);
                }
            }

            while(g < k)
            {
                *g = *l;
                *l = *++g;
                l += cmp(*l, p);
            }
            *g = *l; *l = p;
        }

        // Skip the pivot.
        g = l + (l < high);
        l -= (l > low);

        // Cheaply calculate an
        // eigth of the interval. 
        _8th = x >> 3U;

        // Calculate interval widths.
        
        ls = l - low,
        gs = high - g;

        // If the partition is fairly 
        // balanced, try insertion sort.
        // If insertion sort runtime
        // trends higher than O(n), fall 
        // back to quicksort.
        if(ls >= _8th &&
           gs >= _8th) 
        {
            if(work) goto l1;
            if(!iSort<0,0>(low, l, cmp, low == olow)) 
                goto l1;
            if(!iSort<1,0>(g, high, cmp))
                goto l2;
            goto l3;
        }

        // The partition is not balanced.
        // scramble some elements and
        // try to break the pattern.
        scramble(low, l, ls);
        scramble(g, high, gs);

        // This was a bad partition,
        // so decrement the height.
        // When the height is zero,
        // we will use heapsort.
        --height;

        // Sort left portion.
        l1: if(l > low)
        {
            // If the interval is small
            // enough, use insertion
            // sort.
            if(ls < InsertionThreshold)
                iSort<0,0,0>
                (low, l, cmp, low == olow);

            // If we've seen too many
            // bad partitions, or if
            // the stack overflows,
            // use heapsort.
            else if(height <= 0 || s >= end)
                hSort(low, l, cmp);

            // If the interval isn't sorted
            // push it onto the stack to
            // address later.
            else
            {
                s->a = (void*)low;
                s->b = (void*)l;
                *sh  = uint8_t(height);
                ++s; ++sh;
            }
        }

        // Sort right portion.
        l2: if(high > g) 
        {
            // If the interval is small
            // enough, use insertion
            // sort.
            if(gs < InsertionThreshold)
                iSort<1,0,0>
                (g, high, cmp);

            // If we've seen too many
            // bad partitions, use 
            // heapsort.
            else if(height <= 0)
                hSort(g, high, cmp);

            // If the interval isn't sorted
            // address it now.
            else
            {
                x = gs;
                low = g;
                continue;
            }
        }

        // If there are no intervals
        // left to sort, exit.
        l3: if(s <= stack) return;

        // Pop the next interval
        // from the stack and 
        // address it now.
        --s; --sh;
        low  = (E*) s->a;
        high = (E*) s->b;
        height = int(*sh);
        x = high - low;
    }
}

/**
 * sort 
 * 
 * @tparam E the element type
 * @tparam Cmp the comparator type
 * @param a the array to be sorted
 * @param cnt the size of the the array
 * @param cmp the comparator
 */
template <bool Block = true, typename E, class Cmp>
inline void blipsort
    (
    E* const a,
    const uint32_t cnt,
    const Cmp cmp
    ) 
{
    if(cnt < InsertionThreshold)
    {
        iSort<0,1,0>(a, a + (cnt - 1), cmp);
        return;
    }
    return qSort
    <!std::is_arithmetic<E>::value && 
     !std::is_pointer<E>::value, Block>
        (a, a + (cnt - 1), log2(cnt), cmp);
}}

namespace Arrays 
{
    /**
     * <h1>
     *  <b>
     *  <i>blipsort</i>
     *  </b>
     * </h1> 
     * 
     * <p>
     * Sorts the given array with provided comparator
     * or in ascending order if nonesuch.
     * </p>
     * 
     * @tparam E the element type
     * @tparam Cmp the comparator type
     * @param a the array to be sorted
     * @param cnt the size of the the array
     * @param cmp the comparator
     */
    template <typename E, class Cmp = std::less<>>
    inline void blipsort
        (
        E* const a,
        const uint32_t cnt,
        const Cmp cmp = std::less<>()
        ) 
    {
        Algo::blipsort(a, cnt, cmp);
    }

    /**
     * <h1>
     *  <b>
     *  <i>blipsort_embed</i>
     *  </b>
     * </h1> 
     * 
     * <p>
     * Sorts the given array with provided comparator
     * or in ascending order if nonesuch. Conserves
     * memory by using branchy Hoare instead of 
     * Block Hoare on large types.
     * </p>
     * 
     * @tparam E the element type
     * @tparam Cmp the comparator type
     * @param a the array to be sorted
     * @param cnt the size of the the array
     * @param cmp the comparator
     */
    template <typename E, class Cmp = std::less<>>
    inline void blipsort_embed
        (
        E* const a,
        const uint32_t cnt,
        const Cmp cmp = std::less<>()
        ) 
    {
        Algo::blipsort<0>(a, cnt, cmp);
    }
}

#endif //BMSRS_SORT_H
