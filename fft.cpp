#include <iostream>
#include <complex>
#include <vector>
#include <chrono>
#include <utility>
using namespace std;
using namespace chrono;
//定数は先頭大文字で
const complex<double> I(0,1);
const double Pi = M_PI;

/* TeX format
 * $x^N の根のうち,0 \leq \arg x < \pi を満たす
 * \frac{N}{2}個の根をvectorで返す.$
 */
vector<complex<double> > root_of_unity(int N, bool inverse=false){
	vector<complex<double> > zeta(N);
	zeta[0] = 1;zeta[1] = exp((2*inverse-1)*2*Pi*I/(double)N);
	for(int i=2;i<N/2;i++){
		if(i&1){
			zeta[i] = zeta[1] * zeta[i-1];
		}else{
			zeta[i] = zeta[i>>1] * zeta[i>>1];
		}
	}
	return zeta;
}

//0,1...N-1をビット反転した値をvectorで返す
vector<int> bit_reverse(int N){
	int M = N;
	vector<int> dst(N);
	for(int i=1;i<N;i<<=1){
		M >>= 1;
		for(int j=0;j<i;j++){
			dst[i+j] = dst[j] + M ;
		}
	}
	return dst;
}

/* ifftをするときはinverseフラグを立てる。
 * 値渡しでなく参照渡しにすればもう少し早くなるか?
 */
vector<complex<double> > fft(vector<complex<double> > src, bool inverse=false){
	int N = src.size();
	int m = log2(N);
	vector<complex<double> > dst(N);
	//メモリ節約のためにdstとsrcをポインタでswapする
	vector<complex<double> > *psrc = &src;
	vector<complex<double> > *pdst = &dst;
	vector<complex<double> > zeta = root_of_unity(N, inverse);
	vector<int> rev = bit_reverse(N);
	for(int i=0;i<N;i++){dst[i] = src[rev[i]];}
	for(int j=0;j<m;j++){
		swap(psrc,pdst);
		for(int i=0;i<N;i++){
			//iのj+1ビット目が降りているかどうか
			if(!((i>>j)&1)){
				int k = i|(1<<j);
				//k = iのj+1ビット目を立てた値
				complex<double> tmp = zeta[(i<<(m-1-j))&(1<<(m-1)-1)]*(*psrc)[k];
				(*pdst)[i] = (*psrc)[i] + tmp;
				(*pdst)[k] = (*psrc)[i] - tmp;
			}
		}
	}
	if(!inverse){
		for(int i=0;i<N;i++){(*pdst)[i] = (*pdst)[i]/(double)N;}
	}
	return *pdst;
}

int main(){
	//0,1,..,N-1をfft,ifftして元に戻るか試してみる
	int N = 65536;
	vector<complex<double> > src(N);
	system_clock::time_point  start, end; // 型は auto で可
	for(int i=0;i<N;i++){src[i] = i;}
	start = system_clock::now();
	src = fft(src);
	end = system_clock::now();
	src = fft(src, true);
	for(int i=0;i<N;i++){cout << src[i] << endl;}
	double elapsed = duration_cast<milliseconds>(end-start).count();
	cout << elapsed << "ms" << endl;
}
