// simpix.cpp
// Simulated annealing pixel mapping: rearrange pixels of source image A to match target image B
// using a one-to-one mapping (a permutation).
//
// This version produces BOTH transforms:
//   - A -> B
//   - B -> A
// and makes a 2x2 collage:
//   [1] A     [2] B
//   [3] A->B  [4] B->A

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

/* ---------------- Energy + O(1) swap delta ---------------- */

static long long total_energy(const UInt_t* srcPix, const UInt_t* tgtPix,
                              const vector<int>& p)
{
  long long E = 0;
  const int N = (int)p.size();
  for(int i=0;i<N;i++){
    E += (long long)dist2_rgb(srcPix[p[i]], tgtPix[i]);
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
  double T0     = 4000.0;   // starting temperature
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
  shuffle(p.begin(), p.end(), rng);

  long long E = total_energy(srcPix, tgtPix, p);
  E_init = E;

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

static TASImage build_transform_image(const TASImage& srcImg,
                                      const UInt_t* srcPix,
                                      const vector<int>& perm,
                                      int W, int H)
{
  const int N = W*H;

  TASImage out(srcImg);      // clone => valid internal image
  out.SetEditable(kTRUE);    // allow edits to stick
  UInt_t* outPix = out.GetArgbArray();

  for(int i=0;i<N;i++){
    UInt_t p = srcPix[perm[i]];
    outPix[i] = (p & 0x00FFFFFF) | 0xFF000000; // opaque
  }
  return out;
}

/* ---------------- Main ---------------- */

int main(int argc, char **argv){

  if(argc < 3){
    cout << "Usage: simpix imageA imageB <out_prefix>\n";
    cout << "Outputs: <prefix>_AtoB.png, <prefix>_BtoA.png, collage.png\n";
    return 0;
  }

  TString fA = argv[1];
  TString fB = argv[2];

  TString prefix = (argc > 3) ? argv[3] : "out";
  // allow passing "out.png" too, we'll strip .png if present
  if(prefix.EndsWith(".png")) prefix.ReplaceAll(".png","");

  TString outAtoB = prefix + "_AtoB.png";
  TString outBtoA = prefix + "_BtoA.png";

  cout << "Reading images:\n  A(source)= " << fA << "\n  B(target)= " << fB << "\n";
  cout << "Outputs:\n  " << outAtoB << "\n  " << outBtoA << "\n";

  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  TASImage A(fA.Data());
  TASImage B(fB.Data());

  cout << "ROOT sees:\n";
  cout << "  A: " << A.GetWidth() << " x " << A.GetHeight() << "\n";
  cout << "  B: " << B.GetWidth() << " x " << B.GetHeight() << "\n";

  assert(A.GetWidth()  == B.GetWidth());
  assert(A.GetHeight() == B.GetHeight());

  int W = A.GetWidth();
  int H = A.GetHeight();
  int N = W*H;

  cout << "Pixel Geometry: " << W << " x " << H << "  (N=" << N << ")\n";

  UInt_t* APix = A.GetArgbArray();
  UInt_t* BPix = B.GetArgbArray();

  // Simulated annealing parameters
  SAParams P;
  P.T0 = 4000.0;
  P.alpha = 0.995;
  P.nTemps = 300;
  P.movesPerT_mul = 10;
  P.seed = 12345;

  // ---------- A -> B ----------
  long long E0_AB=0, Ebest_AB=0;
  auto t0 = chrono::high_resolution_clock::now();
  vector<int> permAB = anneal_permutation(APix, BPix, N, P, E0_AB, Ebest_AB);
  auto t1 = chrono::high_resolution_clock::now();
  double secAB = chrono::duration<double>(t1-t0).count();

  cout << "\n[A->B]\n";
  cout << "Initial energy: " << E0_AB << "\n";
  cout << "Best energy:    " << Ebest_AB << "\n";
  cout << "Runtime:        " << secAB << " s\n";

  TASImage imgAtoB = build_transform_image(A, APix, permAB, W, H);
  imgAtoB.WriteImage(outAtoB.Data());

  // ---------- B -> A ----------
  long long E0_BA=0, Ebest_BA=0;
  auto t2 = chrono::high_resolution_clock::now();
  vector<int> permBA = anneal_permutation(BPix, APix, N, P, E0_BA, Ebest_BA);
  auto t3 = chrono::high_resolution_clock::now();
  double secBA = chrono::duration<double>(t3-t2).count();

  cout << "\n[B->A]\n";
  cout << "Initial energy: " << E0_BA << "\n";
  cout << "Best energy:    " << Ebest_BA << "\n";
  cout << "Runtime:        " << secBA << " s\n";

  TASImage imgBtoA = build_transform_image(B, BPix, permBA, W, H);
  imgBtoA.WriteImage(outBtoA.Data());

  // Collage: A, B, A->B, B->A
  TCanvas c("c","simpix",1200,900);
  c.Divide(2,2);
  c.cd(1); A.Draw("X");
  c.cd(2); B.Draw("X");
  c.cd(3); imgAtoB.Draw("X");
  c.cd(4); imgBtoA.Draw("X");
  c.Print("collage.png");

  cout << "\nWrote: " << outAtoB << ", " << outBtoA << ", collage.png\n";
  return 0;
}
