#include "sort.h"

namespace 
{ enum : uint32_t
{
    InsertionThreshold = 88,
    AscendingThreshold = 8,
    LargeDataThreshold = 128
};

/**
 * An interval or subarray.
 */
struct interval final
{ void * a, * b; }; 

/**
 * The DeBruijn constant.
 */
constexpr uint64_t DeBruijn64 =
    0x03F79D71B4CB0A89L;

/**
 * The DeBruijn map from key to integer
 * square index.
 */
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

/**
 * bitScanReverse
 * @authors Kim Walisch
 * @authors Mark Dickinson
 * @param bb bitboard to scan
 * @precondition bb != 0
 * @return index (0..63) of most significant one bit
 */
constexpr int bitScanRev
    (
    uint32_t l
    ) 
{
    assert(l != 0);
    parallelPrefixFill(l);
    return DeBruijnTableF[(int)
        ((l * DeBruijn64) >> 58U)
    ];
}

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
    E const
    el = *i;
    *i = *j;
    *j = el;
}

/**
 *
 * A generic "sift down" method (AKA max-heapify)
 *
 * @tparam E the element type
 * @param a the pointer to the base of the current
 * sub-array
 * @param i the starting index
 * @param size the size of the current sub-array
 */
template<typename E>
inline void siftDown
    (
    E* const a,
    const int i,
    const int size
    ) 
{
    // Store size in
    // a local variable.
    const uint32_t n = size;

    // Establish non-leaf
    // boundary.
    const uint32_t o = n >> 1U;

    // Extract the element
    // to sift.
    E z = a[i];

    // initialize temporary
    // variables.
    uint32_t x = i, l, r;

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
        if(r < n && y < a[r])
            y = a[l = r];

        // if y is less than or
        // equal to the element
        // we are sifting, then
        // we are done.
        if(y <= z) break;

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
 * <b>
 *  <i>Heap Sort</i>
 * </b>
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
template<typename E>
inline void hSort
    (
    E* const low,
    E* const high
    ) 
{
    E* r = high + 1;
    E* const l = low;

    // Build the heap.
    int x = r - l;
    for(int i =
        (x >> 1U); i >= 0; --i)
        siftDown(l, i, x);
    
    // Sort.
    while(l < --r) 
    {
        const E z = *l; *l = *r;
        siftDown(l, 0, --x);
        *r = z;
    }
}

/**
 * <b>
 *  <i>
 * Insertion Sort
 *  </i>
 * </b>
 *
 * <p>
 * Classical ascending insertion sort packaged with a
 * "pairing" optimization to be used in the context of
 * Quick Sort.
 * This optimization is used whenever the portion of
 * the array to be sorted is padded on the left by
 * a portion with lesser elements. The fact that all of
 * the elements on the left are automatically less than
 * the elements in the current portion allows us to skip
 * the costly lower boundary check in the nested loops
 * and insert two elements in one go.
 * </p>
 *
 * @authors Josh Bloch
 * @authors Jon Bently
 * @authors Orson Peters
 * @authors Ellie Moore
 * @tparam E the element type
 * @tparam Are we sorting optimistically?
 * @param leftmost whether this is the leftmost part
 * @param low a pointer to the leftmost index
 * @param high a pointer to the rightmost index
 * left-most partition.
 */
template<typename E, bool Bail = true>
inline bool iSort
    (
    bool leftmost,
    E *const low, 
    E *const high
    ) 
{
    E* l = low;
    E* r = high;
    int moves = 0;
    if (leftmost) 
    {
        // Traditional
        // insertion
        // sort.
        for (E *i = l + 1; i <= r; ++i) 
        {
            E t = *i, *j = i - 1;
            for (; j >= l && t < *j; --j)
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
        // Pair insertion sort.
        // Skip elements that are
        // in ascending order.
        do if (l++ >= r) return true;
        while (*l >= *(l - 1));

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
            if (ey < ex) 
            {
                ex = ey;
                ey = *i;
                ++moves;
            }

            // Insert the two
            // in one downward
            // motion.
            while (ey < *--i)
                i[2] = *i;
            (++i)[1] = ey;
            while (ex < *--i)
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
        while (ez < *--r)
            r[1] = *r;
        r[1] = ez;
    }
    return true;
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
inline void scramble
    (
    E* const low, 
    E* const high, 
    const uint32_t len
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
 * The sort stack(s). 
 */
interval stack[bitScanRev(UINT32_MAX)];
uint8_t heights[bitScanRev(UINT32_MAX)];
interval * end = stack + bitScanRev(UINT32_MAX);

/**
 * <b>
 *  <i>Blipsort</i>
 * </b>
 *
 * <p>
 * See readme (too lazy to type this right now)
 * </p>
 *
 * @authors Josh Bloch
 * @authors Jon Bently
 * @authors Orson Peters
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
<typename E, bool Root = true>
void qSort
    (
    E * low,
    E * high,
    int height
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
        if(*low <= *cl  |
           *cl  <= *sl  |
           *sl  <= *mid |
           *mid <= *sr  |
           *sr  <= *cr  |
           *cr  <= *high)
        {
            
            if(*low  < *cl) 
                cl = low;
            if(*high > *cr) 
                cr = high;

            if (*sl < *cl) 
            {
                E e = *sl;
                *sl = *cl;
                *cl =   e;
            }

            if (*mid < *sl) 
            {
                E e  = *mid;
                *mid =  *sl;
                *sl  =    e;
                if (e < *cl) 
                {
                    *sl = *cl;
                    *cl =   e;
                }
            }

            if (*sr < *mid) 
            {
                E e  =  *sr;
                *sr  = *mid;
                *mid =    e;
                if (e < *sl) 
                {
                    *mid = *sl;
                    *sl  =   e;
                    if (e < *cl) 
                    {
                        *sl = *cl;
                        *cl =   e;
                    }
                }
            }

            if (*cr < *sr) 
            {
                E e = *cr;
                *cr = *sr;
                *sr =   e;
                if (e < *mid) 
                {
                    *sr  = *mid;
                    *mid =    e;
                    if (e < *sl) 
                    {
                        *mid = *sl;
                        *sl  =   e;
                        if (e < *cl) 
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
        // is no big deal for
        // branchless Lomuto.
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
            if(h == *sl  || 
               h == *mid || 
               h == *sr) 
            {
                l = low - 1;
                g = high + 1;

                // skip over data
                // in place.         
                while(*--g > h);

                E e = *g;
                *g = h + 1;
                while(*++l == h);
                *g = e;
                
        /**
         * Partition left by branchless Lomuto scheme
         * 
         * During partitioning:
         * 
         * +-------------------------------------------------------------+
         * |  ... == p  |  ... > p  | * |     ... ? ...      |  ... > p  |
         * +-------------------------------------------------------------+
         * ^            ^           ^                        ^           ^
         * low          l           k                        g         high
         * 
         * After partitioning:
         * 
         * +-------------------------------------------------------------+
         * |           ... == p           |            > p ...           |
         * +-------------------------------------------------------------+
         * ^                              ^                              ^
         * low                            l                           high
         */
                k = l; p = *l;
                while(k < g)
                {
                    *k++ = *l;
                    *l = *k;
                    l += (*l == h);
                }
                *k = *l;
                *l = p;
                l += (p == h);
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
                    iSort<E, false>
                    (false, low, high);
                    goto l3;
                }

                // Continue with quicksort.
                continue;
            }
        }

        // Initialize l and g.
        // "less" and "great"
        // respectively.
        l = low - 1, k = high + 1;

        // Assign midpoint to pivot
        // variable.
        p = *mid;

        // skip over data
        // in place.
        while(*++l < p);

        // Bring left end inside.
        // Left end will be
        // replaced and pivot will
        // be swapped back later.
        *mid = *l;

        // Avoid running past low
        // end. place a stopper in
        // the gap.
        *l = p - 1;

        // skip over data
        // in place.
        while(*--k >= p);

        // Will we do a significant 
        // amount of work during 
        // partitioning?
        work = 
        ((l - low) + (high - k)) 
            < (x >> 1U);

        g = l;

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
        while(g < k)
        {
            *g++ = *l;
            *l = *g;
            l += (*l < p);
        }
        *g = *l; *l = p;

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
            if(!iSort(low == olow, low, l)) 
                goto l1;
            if(!iSort(false, g, high))
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
            // If we've seen too many
            // bad partitions, or if
            // the stack overflows,
            // use heapsort.
            if(height <= 0 || s >= end)
                hSort(low, l);

            // If the interval is small
            // enough, use insertion
            // sort.
            else if(ls < InsertionThreshold)
                iSort<E, false>
                (low == olow, low, l);

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
            // If we've seen too many
            // bad partitions, use 
            // heapsort.
            if(height <= 0)
                hSort(g, high);

            // If the interval is small
            // enough, use insertion
            // sort.
            else if(gs < InsertionThreshold)
                iSort<E, false>
                (false, g, high);

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
}}

namespace Arrays 
{
    template <typename E>
    inline void blipsort
        (
        E* const a,
        const uint32_t cnt
        ) 
    {
        if(cnt < InsertionThreshold)
        {
            iSort<E, false>(true, a, a + (cnt - 1));
            return;
        }

        // floor of log base 2 of cnt.
        const int log2Cnt = bitScanRev(cnt);
        return qSort(a, a + (cnt - 1), log2Cnt);
    }

    template void
    blipsort<int64_t>(int64_t*, uint32_t);
    template void
    blipsort<int32_t>(int32_t*, uint32_t);
    template void
    blipsort<int16_t>(int16_t*, uint32_t);
    template void
    blipsort<int8_t>(int8_t*, uint32_t);
}

#include <algorithm>
#include <iostream>
#include <chrono>
#include <vector>
#include "pdqsort.h"
#ifdef _WIN32
    #include <intrin.h>
    #define rdtsc __rdtsc
#else
    #ifdef __i586__
        static __inline__ unsigned long long rdtsc() {
            unsigned long long int x;
            __asm__ volatile(".byte 0x0f, 0x31" : "=A" (x));
            return x;
        }
    #elif defined(__x86_64__)
        static __inline__ unsigned long long rdtsc(){
            unsigned hi, lo;
            __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
            return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
        }
    #else
        #error no rdtsc implementation
    #endif
#endif

void shellsort(int64_t* a, int l, int r)
{ int i, j, h; int64_t v;
int incs[16] = { 1391376, 463792, 198768, 86961, 33936, 13776, 4592, 1968, 861, 336, 112, 48, 21, 7, 3, 1 };
for ( int k = 0; k < 16; k++) {
    for (h = incs[k], i = l+h; i <= r; i++) {
        v = a[i]; j = i;
        while (j >= h && a[j-h] > v)
        { a[j] = a[j-h]; j -= h; }
        a[j] = v; 
    }
}
}

struct pair { const char* name; int64_t* (* func)(int64_t* t, int size); };

void shuffle(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i)
    {
        const int r = rand() % size;
        const uint64_t temp = t[r];
        t[r] = t[i];
        t[i] = temp;
    }
} 

int64_t* shuffled_many_equal(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i) 
        t[i] = rand() % (size/100);
    return t;
}

int64_t* shuffled_some_equal(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i) 
        t[i] = rand() % (size/10);
    return t;
}

int64_t* shuffled(int64_t* t, int size) 
{
    for (int i = 0; i < size; ++i) t[i] = i;
    shuffle(t, size);
    return t;
}

int64_t* shuffled16(int64_t* t, int size) 
{
    for (int i = 0; i < size; ++i) t[i] = (i + 1) % 16;
    shuffle(t, size);
    return t;
}

int64_t* all_equal(int64_t* t, int size)  
{
    for (int i = 0; i < size; ++i) t[i] = 0;
    return t;
}

int64_t* ascending(int64_t* t, int size)   
{
    for (int i = 0; i < size; ++i) t[i] = i;
    return t;
}

int64_t* descending(int64_t* t, int size)   
{
    for (int i = 0; i < size; ++i) t[i] = size - i;
    return t;
}

int64_t* ascending_saw(int64_t* t, int size)   
{
    for (int i = 0; i < size; ++i) t[i] = (i + 1) % (size / 25);
    return t;
}

int64_t* descending_saw(int64_t* t, int size)   
{
    for (int i = 0; i < size; ++i) t[i] = (size - i) % (size / 25);
    return t;
}

int64_t* pipe_organ(int64_t* t, int size)   
{
    for (int i = 0; i < size/2; ++i) t[i] = i;
    for (int i = size/2; i < size; ++i) t[i] = size - i;
    return t;
}

int64_t* push_front(int64_t* t, int size)
{
    int i = 1;
    for (; i < size; ++i) t[i - 1] = i;
    t[i - 1] = 0;
    return t;
}

int64_t* push_middle(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i) {
        if (i != size/2) t[i] = i;
        else t[i] = 0;
    }
    return t;
}

int64_t* descending_blinds(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i) {
        if (i % 2 != 0) t[i] = size - i;
        else t[i] = 0;
    }
    return t;
}

int64_t* descending_madness(int64_t* t, int size)
{
    for (int i = 0, j = 20; i < size; ++i) {
        t[i] = (size - i);
        if(j++ % 200 != 0) --i;
    }
    return t;
}

int64_t* descending_and_ascending(int64_t* t, int size)
{
    for (int i = 0; i < size; ++i) {
        t[i] = i;
    }

    for (int i = 0; i < size; i += 2) {
        t[i] = size-i;
    }
    return t;
}

int64_t* blipsort_adaptor(int64_t* t, int size)
{
    Arrays::blipsort<int64_t>(t, size);
    return nullptr;
}

int64_t* pdqsort_adaptor(int64_t* t, int size)
{
    pdqsort<int64_t*>(t, t + size);
    return nullptr;
}

int64_t* stdsort_adaptor(int64_t* t, int size)
{
    std::sort(t, t + (size - 1));
    return nullptr;
}

int64_t* shellsort_adaptor(int64_t* t, int size)
{
    shellsort(t, 0, size - 1);
    return nullptr;
}

int64_t* hsort_adaptor(int64_t* t, int size)
{
    hSort(t, t + (size - 1));
    return nullptr;
}

pair s[] = 
{
    {"shuffled", shuffled},
    {"shuffled16", shuffled16},
    {"shuffled_many_equal", shuffled_many_equal},
    {"shuffled_some_equal", shuffled_some_equal},
    {"all_equal", all_equal},
    {"ascending", ascending},
    {"descending", descending},
    {"pipe_organ", pipe_organ},
    {"push_front", push_front},
    {"push_middle", push_middle},
    // {"descending_blinds", descending_blinds},
    // {"descending_madness", descending_madness},
    // {"ascending_saw", ascending_saw},
    // // {"descending_saw", descending_saw}
    {"descending and ascending", descending_and_ascending}
};

pair sorts[] = 
{
    {"blipsort", blipsort_adaptor},
    {"stdsort", stdsort_adaptor},
    {"pdqsort", pdqsort_adaptor},
    // {"shellsort", shellsort_adaptor},
    // {"heapsort", hsort_adaptor}
};

int main()
{
    // int64_t arr[500];
    // for(int i = 0; i < 500; ++i) arr[i] = 500 - i;
    
    // Arrays::uSort<int64_t>(arr, 500);
    // for(int i = 0; i < 500; ++i) std::cout << arr[i] << ' ';
    // std::cout << '\n';

    // return 0;

    for(pair& dist: s)
    {
    // pair dist = s[6];
        srand(time(NULL));
        for(pair& sort_ : sorts )
        {
            for(int size: {1000000})
            {
                int64_t * t = new int64_t[size];
                std::chrono::time_point<std::chrono::high_resolution_clock> total_start, total_end;
                std::vector<uint64_t> cycles;

                total_start = std::chrono::high_resolution_clock::now();
                total_end = std::chrono::high_resolution_clock::now();
                while (std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start).count() < 5000) {
                    t = dist.func(t, size);
                    uint64_t start = rdtsc();
                    sort_.func(t, size);
                    uint64_t end = rdtsc();
                    cycles.push_back(uint64_t(double(end - start) / size + 0.5));
                    total_end = std::chrono::high_resolution_clock::now();
                    // if (!std::is_sorted(v.begin(), v.end())) {
                    //     std::cerr << "sort failed: ";
                    //     std::cerr << size << " " << distribution.first << " " << sort.first << "\n";
                    // }
                }

                std::sort(cycles.begin(), cycles.end());

                std::cerr << size << " " << dist.name << " " << sort_.name
                          << " " << cycles[cycles.size()/2] << "\n";
                std::cout << size << " " << dist.name << " " << sort_.name
                          << " " << cycles[cycles.size()/2] << "\n";  
                delete[] t;
            }
        }
    }
}

// #include<iostream>
// #include<ctime>
// #include<stdio.h>
// #include<chrono>
// #include<algorithm>
// #include <bitset>
// #include <climits>
// #include <cstring>
// #include <iostream>
// #include "pdqsort.h"

// // int main() {
// //     int a[1000000];
// //     srand(time(NULL));
// //     for(int i = 0; i < 1000000; ++i)
// //         a[i] = rand() % 1000000;
// //     Arrays::uSort<int32_t>(a, 1000000);
// //     for(int i = 1; i < 1000000; ++i)
// //     {
// //         if( a[i - 1] > a[i] )
// //         {
// //             std::cout << "\nNOT SORTED\n";
// //             return 0;
// //         }
// //     }

// //     std::cout << "\nSORTED\n";
// //     return 0;       
// // }


// /**
//  * A stack of intervals, to
//  * simulate recursion. 
//  */

// constexpr int num = 100000;
// int b[10000000];
// int main()
// {
//     srand(time(NULL));
//     for(int o = 0; o < 100; ++o)
//     {
//     typedef std::chrono::high_resolution_clock Clock;
//     uint64_t t = 0;
//     int a[num];
//     for(int c = 0; c < 500; ++c)
//     {
//         std::memset(b, 0, sizeof(b));
//         for(int i = 0; i < num; ++i) {
//             b[a[i] = rand() % num]++;
//             // std::cout << a[i] << ' ';
//         }
//         // std::cout << '\n';

//         // for(int i = 0; i < num; ++i) {
//         //     if(i % 1000)
//         //         a[i] = rand() % num;
//         //     else a[i] = i;
//         // }
//         auto t1 = Clock::now();
//         Arrays::blipsort<int32_t>(a, num);
//         // std::sort(a, a + num);
//         auto t2 = Clock::now();
//         t += std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
//         int i;
//         for(i = 1; i < num; ++i)
//         {
//             // std::cout << a[i - 1] << ' ';
//             if(a[i - 1] > a[i]) 
//             {
//                 std::cout << "NOT SORTED!!!";
//                 goto out;
//             }
//         }
//         // std::cout << a[i - 1] << '\n';

//         for(int i = 0; i < num; ++i)
//         {
//             b[a[i]]--;

//         }

//         for(int j = 0; j < num; ++j)
//         {
//             if(b[a[j]] != 0)
//             {
//                 std::cout << "MISSING ELEMENTS!!!\n";
//                 std::cout << b[a[j]]  << '\n';
//                 goto out;
//             }
//         }
//     }
//     std::cout << t << '\n';
//     out: continue;
//     }
//     return 0;
// }

// // int main()
// // {
// //     int a[] = {0,0,1,2,3,4,5,9,10,11};

// //     int* u = a;
// //     int* q = a + 9;
// //     int* mid = a + (((a + 9) - a) >> 1U);
// //     while(u < mid)
// //     {
// //         int e = *u;
// //         *u++ = *q;
// //         *q-- = e;
// //     }

// //     // if(Arrays::iSort<int, true>(a + 0, a + 7))
// //     //     std::cout << "true\n";
// //     // else std::cout << "false\n";
// //     for(int i = 0; i < 10; ++i)
// //         std::cout << a[i] << ' ';
// // }