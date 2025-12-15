// simpix.cpp
// Simulated annealing pixel mapping: rearrange pixels of source image A to match target image B
// using a one-to-one mapping (a permutation).

#include "TROOT.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"

#include <cassert>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

using namespace std;

/* ---------------- Pixel utilities ---------------- */

static inline int R(UInt_t p){ return (p >> 16) & 0xFF; }
static inline int G(UInt_t p){ return (p >>  8) & 0xFF; }
static inline int B(UInt_t p){ return (p      ) & 0xFF; }

static inline int dist2_rgb(UInt_t a, UInt_t b){
  int dr = R(a) - R(b);
  int dg = G(a) - G(b);
  int db = B(a) - B(b);
  return dr*dr + dg*dg + db*db;
}

// reverse entire buffer (same flip you had originally) => 180-degree rotation
static void flip_reverse_buffer(UInt_t* pix, long long N){
  for(long long i=0;i<N/2;i++){
    UInt_t tmp = pix[i];
    pix[i] = pix[N - i - 1];
    pix[N - i - 1] = tmp;
  }
}

/* ---------------- Energy + O(1) swap delta ---------------- */

static long long total_energy(const UInt_t* srcPix, const UInt_t* tgtPix,
                              const vector<int>& p)
{
  long long E = 0;
  const int N = (int)p.size();
  for(int i=0;i<N;i++){
    E += dist2_rgb(srcPix[p[i]], tgtPix[i]);
  }
  return E;
}

// Delta-E for swapping p[i] and p[j]
static inline long long delta_swap(const UInt_t* srcPix, const UInt_t* tgtPix,
                                  const vector<int>& p, int i, int j)
{
  const int pi = p[i];
  const int pj = p[j];

  long long before =
      (long long)dist2_rgb(srcPix[pi], tgtPix[i]) +
      (long long)dist2_rgb(srcPix[pj], tgtPix[j]);

  long long after  =
      (long long)dist2_rgb(srcPix[pj], tgtPix[i]) +
      (long long)dist2_rgb(srcPix[pi], tgtPix[j]);

  return after - before;
}

/* ---------------- Simulated annealing ---------------- */

struct SAParams {
  double T0     = 4000.0;   // starting temperature (tune if needed)
  double alpha  = 0.995;    // cooling factor per temperature step
  int    nTemps = 300;      // number of temperature steps
  int    movesPerT_mul = 10; // movesPerT = movesPerT_mul * N
  unsigned seed = 12345;
};

static vector<int> anneal_permutation(const UInt_t* srcPix, const UInt_t* tgtPix,
                                      int N, const SAParams& P,
                                      long long& E_init, long long& E_best)
{
  mt19937 rng(P.seed);
  uniform_int_distribution<int> pick(0, N-1);
  uniform_real_distribution<double> uni(0.0, 1.0);

  // permutation p: output position i gets source pixel index p[i]
  vector<int> p(N);
  for(int i=0;i<N;i++) p[i]=i;

  // random initial mapping (better exploration)
  shuffle(p.begin(), p.end(), rng);

  long long E = total_energy(srcPix, tgtPix, p);
  E_init = E;

  // keep track of best seen permutation
  vector<int> p_best = p;
  E_best = E;

  double T = P.T0;
  const int movesPerT = P.movesPerT_mul * N;

  for(int t=0; t<P.nTemps; t++){
    int accepted = 0;

    for(int m=0; m<movesPerT; m++){
      int i = pick(rng);
      int j = pick(rng);
      if(i==j) continue;

      long long dE = delta_swap(srcPix, tgtPix, p, i, j);

      if(dE <= 0){
        swap(p[i], p[j]);
        E += dE;
        accepted++;
      } else {
        double prob = exp(-(double)dE / T);
        if(uni(rng) < prob){
          swap(p[i], p[j]);
          E += dE;
          accepted++;
        }
      }

      if(E < E_best){
        E_best = E;
        p_best = p;
      }
    }

    if(t % 25 == 0){
      double acc = (double)accepted / (double)movesPerT;
      cout << "step " << t
           << "  T=" << T
           << "  E=" << E
           << "  best=" << E_best
           << "  acc=" << acc << "\n";
    }

    T *= P.alpha;
  }

  return p_best;
}

/* ---------------- Main ---------------- */

int main(int argc, char **argv){

  if(argc < 3){
    cout << "Usage: simpix imageA imageB <out_AtoB.png>\n";
    return 0;
  }

  TString fsrc = argv[1];
  TString ftgt = argv[2];
  TString fout = (argc > 3) ? argv[3] : "out_AtoB.png";

  TString fout_flipped = fout;
  if(fout_flipped.EndsWith(".png")) fout_flipped.ReplaceAll(".png","_flipped.png");
  else fout_flipped += "_flipped.png";

  cout << "Reading images:\n  A(source)= " << fsrc << "\n  B(target)= " << ftgt << "\n";
  cout << "Output: " << fout << "\n";

  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  TASImage src(fsrc.Data());
  TASImage tgt(ftgt.Data());

  cout << "ROOT sees:\n";
  cout << "  src: " << src.GetWidth() << " x " << src.GetHeight() << "\n";
  cout << "  tgt: " << tgt.GetWidth() << " x " << tgt.GetHeight() << "\n";

  assert(src.GetWidth()  == tgt.GetWidth());
  assert(src.GetHeight() == tgt.GetHeight());

  int W = src.GetWidth();
  int H = src.GetHeight();
  int N = W*H;

  cout << "Pixel Geometry: " << W << " x " << H << "  (N=" << N << ")\n";

  UInt_t* srcPix = src.GetArgbArray();
  UInt_t* tgtPix = tgt.GetArgbArray();

  // Simulated annealing parameters (tune if needed)
  SAParams P;
  P.T0 = 4000.0;
  P.alpha = 0.995;
  P.nTemps = 300;
  P.movesPerT_mul = 10;
  P.seed = 12345;

  long long E0=0, Ebest=0;

  auto t0 = chrono::high_resolution_clock::now();
  vector<int> bestPerm = anneal_permutation(srcPix, tgtPix, N, P, E0, Ebest);
  auto t1 = chrono::high_resolution_clock::now();

  double seconds = chrono::duration<double>(t1-t0).count();
  cout << "Initial energy: " << E0 << "\n";
  cout << "Best energy:    " << Ebest << "\n";
  cout << "Runtime:        " << seconds << " s\n";

  // Build output image from best permutation
  TASImage out(W, H);
  UInt_t* outPix = out.GetArgbArray();

  for(int i=0;i<N;i++){
    UInt_t p = srcPix[bestPerm[i]];
    outPix[i] = (p & 0x00FFFFFF) | 0xFF000000; // force opaque
  }

  // Save correct orientation output
  out.WriteImage(fout.Data());

  // Make flipped version (reverse-buffer flip = 180-degree rotation)
  // TASImage outFlip(out);
  // UInt_t* outFlipPix = outFlip.GetArgbArray();
  // flip_reverse_buffer(outFlipPix, (long long)N);
  // outFlip.WriteImage(fout_flipped.Data());

  TASImage outFlip(out);
  outFlip.Flip(180);                 // 180-degree rotation
  outFlip.WriteImage(fout_flipped.Data());

  // Collage: A, B, out, outFlipped
  TCanvas c("c","simpix",1200,900);
  c.Divide(2,2);
  c.cd(1); src.Draw("X");
  c.cd(2); tgt.Draw("X");
  c.cd(3); out.Draw("X");
  c.cd(4); outFlip.Draw("X");
  c.Print("collage.png");

  cout << "Wrote: " << fout << ", " << fout_flipped << ", collage.png\n";
  return 0;
}
