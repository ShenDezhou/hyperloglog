#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

class HyperLogLog {
  public:
    // 4 <= b <= 16
    explicit HyperLogLog(int64_t b) : b_(b), m_(1ll << b), ms_(1ll << b) {}

    void update(uint64_t value) {
        int64_t j = value & ((1ll << b_) - 1ll);
        uint64_t zeros = nlz5(value >> b_) - b_ + 1ll;
        ms_[j] = std::max(ms_[j], zeros);
    }

    double estimate() const {
        double e = raw_estimate();

        if (e <= 2.5*m_) {
            int zeros = num_zeros();
            if (zeros != 0) {
                e = m_*log((double)m_/zeros);
            }
        } else if (e > pow(2ll, 32)/30.0) {
            e = log1p(e*-1/pow(2ll, 32))*pow(2ll, 32)*-1ll;
        }
        return e;
    }

    double raw_estimate() const {
        double e = 0;
        for(std::vector<uint64_t>::const_iterator iter=ms_.begin();iter!=ms_.end();++iter) {
            e += 1.0/(1ll << *iter);
        }
        e = 1.0/e;
        e *= alpha()*m_*m_;

        return e;
    }

    double alpha() const {
        switch (m_) {
            case 16:
                return 0.673;
            case 32:
                return 0.697;
            case 64:
                return 0.709;
            default:
                return 0.7213/(1ll + 1.079/m_);
        }
    }

    int num_zeros() const {
        int count = 0;
        for(std::vector<uint64_t>::const_iterator iter=ms_.begin();iter!=ms_.end();++iter) {
            if (*iter == 0)
                count++;
        }
        return count;
    }

    std::map<int64_t, int64_t>* histogram() const {
        std::map<int64_t, int64_t>* h = new std::map<int64_t, int64_t>;

        for(std::vector<uint64_t>::const_iterator iter=ms_.begin();iter!=ms_.end();++iter) {
            if (h->find(*iter) != h->end()) {
                h->at(*iter)++;
            } else {
                (*h)[*iter] = 1ll;
            }
        }

        return h;
    }

    int64_t b() const { return b_; }

  private:
    // These from the book "Hacker's Delight"
    uint64_t pop(uint64_t x) {
        x = x - ((x >> 1) & 0x5555555555555555ll);
        x = (x & 0x3333333333333333ll) + ((x >> 2) & 0x3333333333333333ll);
        x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fll;
        x = x + (x << 8);
        x = x + (x << 16);
        return x >> 24;
    }

    uint64_t nlz5(uint64_t x) {
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >>16);
        return pop(~x);
    }

    std::vector<uint64_t> ms_;
    int64_t b_;
    int64_t m_;
};
