#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

int main(void) {
    vector<int> myvector;
    for(int i = 0; i < 5; ++i){
        myvector.push_back(i);
    }

    random_shuffle(myvector.begin(), myvector.end()); 
    for (vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it) {
        cout << ' ' << *it;
    }

    cout << '\n';

    random_shuffle(myvector.begin(), myvector.end()); 
    for (vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it) {
        cout << ' ' << *it;
    }
    return 1;
}

