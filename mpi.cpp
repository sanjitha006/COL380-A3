#include <mpi.h>
#include <bits/stdc++.h>
#define ll long long
#define pb push_back
using namespace std;

enum Msgtags{
    GIVE_TASK=1,
    REQ_TASK=2,
    TERMINATE=3,
    UPDATE_P=4,
    RESULT=5,
}

static int N,E;
static ll B;
static vector<int>profit;
static vector<int>cost;
static vector<uint8_t>adj;
inline bool edge(int u,int v){
    return adj[u*n+v];
}

static vector<int> pc_sort;
static int pmax=0;
static vector<int>best;

struct Problem{
    vector<int>cand;
    vector<int>chosen;
    ll p=0;
    ll w=0;
}

static vector<int> makemsg(const Problem& s){
    vector<int>msg;
    b.pb(s.cand.size());
    b.insert(b.end(),s.cand.begin(),s.cand.end());
    b.pb(s.chosen.size());
    b.insert(b.end(),s.chosen.begin(),s.chosen.end());
    b.pb(s.p);
    b.pb(s.w);
    return b;
}

static Problem unpackmsg(const int* d){
    Problem s;
    int i=0;
    int cs=d[i++];
    s.cand.assign(d+i,d+i+cs);
    i+=cs;
    int ch=d[i++];
    s.chosen.assign(d+i,d+i+ch);
    i+=ch;
    i++;
    s.p=d[i];
    i++;
    s.w=d[i];
    return s;
}

static int structbound(const vector<int>& cand){
    int n=cand.size();
    if(n==0) return 0;
    vector<int>colour(n,-1);
    int ncolors=0;
    int found=0;
    int f=1;
    int res=0;
    for(int i=0;i<n;i++){
        found=0;
        for(int c=0;c<ncolors;c++){
            f=1;
            for(int j=0;j<i;j++){
                if(colour[j]==c and edge(cand[i],cand[j])){
                    f=0;break;
                }
            }
            if(f==0) continue;
            else{
                colout[i]=c;
                found=1;
                break;
            }
        }
        if(found==0){
            ncolors++;
            res+=profit[cand[i]];
            //is cand sorted in order? in terms of the profit...
            // //then only directly add profit of new member
        }
    }
    return res;

}

static int knapbound(const vector<int>&cand,int budget){
    if(budget<=0 or cand.empty()) return 0;
    vector<bool>pres(N,false);
    for(int v:cand) pres[v]=true;
    int res=0;
    int currb=0;
    int total=cand.size();
    int seen=0;
    for(int v:pc_sort){
        if(seen==total) break;
        if(pres[v]){
            seen++;
            if((currb+cost[v])<=budget){
                res+=profit[v];
                currb+=cost[v];
            }
            else{
                res+=((profit[v]*(budget-currb))/cost[v]);
                return res;
            }
        }
    }
    return res;
}

static void branch(const Problem& s,stack<Problem>&st){
    //prune
    if((s.p+structbound(s.cand))<=pmax) return;
    if((s,p+knapbound(s.cand,B-s.W))<=pmax) return;

    for(int i=cand.size()-1;i>=0;i--){
        s.cand[i];
        if((s.w+cost[s.cand[i]])<=B){
            if((s.p+profit[s.cand[i]])>pmax){
                pmax=s.p+profit[s.cand[i]];
                best=s.chosen;
                best.pb(s.cand[i]);
            }
            Problem child;
            child.p=s.p+profit[s.cand[i]];
            child.w=s.w+cost[s.cand[i]];
            child.chosen=s.chosen;
            child.chosen.pb(s.cand[i]);
            child.cand={};
            for(int j=i+1;j<cand.size();j++){
                if(edge(s.cand[i],s.cand[j])){
                    child.cand.pb(s.cand[j]);
                }
            }
            st.pb(move(child));
        }
    }
}

