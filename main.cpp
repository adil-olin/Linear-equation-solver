#include <bits/stdc++.h>
using namespace std;

#ifdef OL_IN_ONE
#include "debug.h"
#else
#define debug(args...)
#endif

#define eps 0.000001

struct MATRIX
{
    int n,m;
    bool exist;
    vector<vector<long double>>val;
    set<string>variable;
    map<string,string>ans;

    MATRIX(int _n,int _m)
    {
        n=_n;
        m=_m;
        val.resize(n,vector<long double>(m,0));
        ans.clear();
        exist = true;
    }
    MATRIX(int _n,int _m,set<string>var)
    {
        n=_n;
        m=_m;
        val.resize(n,vector<long double>(m,0));
        variable=var;
        ans.clear();
        exist=true;
    }
    void print()
    {
        cout<<"[";
        int xx=0;
        for(auto x:variable)
        {
            if(xx!=0)
            {
                cout<<" ";
            }
            cout<<x;
            xx++;
        }
        cout<<"]"<<endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                cout << val[i][j] << " ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
};

string make_normal(string s)
{
    string ans;
    for (int i = 0; i < (int)s.size(); i++)
    {
        if (s[i] == '+' || s[i] == '-' || s[i] == '=')
        {
            ans = ans + ' ' + s[i] + ' ';
        }
        else if (i + 1 < (int)s.size() && s[i] <= '9' && s[i] >= '0' && s[i + 1] != '.' && !(s[i+1]>='0' && s[i+1]<='9'))
        {
            ans = ans + s[i] + ' ';
        }
        else
        {
            ans = ans + s[i];
        }
    }
    return ans;
}

string make_parameter(int n)
{
    int x=n%11,y=n/11; // p,q,r,s.....z (there are 11 of them)
    string ans,tmp;
    char c='p'+x;
    ans+=c;
    if(y)
    {
        tmp=to_string(y);
    }
    ans=ans+tmp;
    return ans;
}

MATRIX make_row_form(MATRIX mat)
{
    vector<pair<int,int>>row(mat.n);
    for(int i=0;i<mat.n;i++)
    {
        int cnt=0;
        for(int j=0;j<mat.m;j++)
        {
            if(abs(mat.val[i][j])<eps)
            {
                cnt++;
            }
            else
            {
                row[i].first=cnt;
                row[i].second=i;
                break;
            }
        }
    }
    sort(row.begin(),row.end());
    MATRIX ans(mat.n,mat.m,mat.variable);
    for(int i=0;i<mat.n;i++)
    {
        for(int j=0;j<mat.m;j++)
        {
            ans.val[i][j]=mat.val[row[i].second][j];
        }
    }
    return ans;
}

MATRIX make_row_echelon(MATRIX mat)
{
    mat=make_row_form(mat);
    for(int i=0;i<mat.n;i++)
    {
        int id=-1;
        for(int j=0;j<mat.m;j++)
        {
            if(abs(mat.val[i][j])>eps)
            {
                id=j;
                break;
            }
        }
        if(id==-1)
        {
            continue;
        }
        for(int ii=i+1;ii<mat.n;ii++)
        {
            long double tmp=mat.val[ii][id]/mat.val[i][id];
            for(int j=0;j<mat.m;j++)
            {
                mat.val[ii][j]-=(mat.val[i][j]*tmp);
                if(abs(mat.val[ii][j])<eps)
                {
                    mat.val[ii][j]=0;
                }
            }
        }
    }
    for(int i=0;i<mat.n;i++)
    {
        long double id=1;
        for(int j=0;j<mat.m;j++)
        {
            if(abs(mat.val[i][j])>eps)
            {
                id=mat.val[i][j];
                break;
            }
        }
        for(int j=0;j<mat.m;j++)
        {
            if(abs(mat.val[i][j])>eps)
            {
                mat.val[i][j]/=id;
            }
        }
    }
    return mat;
}

MATRIX make_reduced_row_echelon(MATRIX mat)
{
    for(int i=0;i<mat.n;i++)
    {
        int id=-1;
        for(int j=0;j<mat.m-1;j++)
        {
            if(abs(mat.val[i][j])>eps)
            {
                id=j;
                break;
            }
        }
        if(id==-1)
        {
            continue;
        }
        for(int ii=0;ii<mat.n;ii++)
        {
            if(ii==i)
            {
                continue;
            }
            if(abs(mat.val[ii][id])>eps)
            {
                long double tmp=mat.val[ii][id]/mat.val[i][id];
                for(int j=0;j<mat.m;j++)
                {
                    mat.val[ii][j]-=(tmp*mat.val[i][j]);
                    if(abs(mat.val[ii][j])<eps)
                    {
                        mat.val[ii][j]=0;
                    }
                }
            }
        }
    }
    return mat;
}

MATRIX make_ans(MATRIX mat)
{
    mat.ans.clear();
    int x=0;
    string got[mat.m];
    long double chk[mat.m];
    memset(chk,0,sizeof(chk));
    for(int i=0;i<mat.n;i++)
    {
        int cnt=0;
        for(int j=0;j<mat.m;j++)
        {
            if(abs(mat.val[i][j])<eps)
            {
                cnt++;
            }
        }
        if(cnt==mat.m-1 && abs(mat.val[i][mat.m-1])>eps)
        {
            mat.exist=false;
        }
    }
    for(int i=mat.n-1;i>=0;i--)
    {
        int id=-1;
        for(int j=0;j<mat.m;j++)
        {
            if(mat.val[i][j]!=0)
            {
                id=j;
                break;
            }
        }
        if(id==-1)
        {
            continue;
        }
        chk[id]=mat.val[i][mat.m-1];
        string tmp_ans;
        for(int j=mat.m-2;j>id;j--)
        {
            chk[id]-=(mat.val[i][j]*chk[j]);
            if(mat.val[i][j]==0)
            {
                continue;
            }
            if(mat.val[i][j]>0)
            {
                tmp_ans+=" - ";
                if(!got[j].empty())
                {
                    string tmp2=to_string(mat.val[i][j]);
                    tmp_ans+=tmp2;
                    tmp_ans+=got[j];
                    continue;
                }
                while (mat.variable.find(make_parameter(x))!=mat.variable.end())
                {
                    x++;
                }
                string tmp=make_parameter(x),tmp2;
                x++;
                tmp2=to_string(mat.val[i][j]);
                got[j]=tmp;
                if(abs(mat.val[i][j])!=1)
                {
                    tmp_ans+=tmp2;
                }
                tmp_ans+=tmp;
            }
            else
            {
                tmp_ans+=" + ";
                mat.val[i][j]=-mat.val[i][j];
                if(!got[j].empty())
                {
                    string tmp2=to_string(mat.val[i][j]);
                    tmp_ans+=tmp2;
                    tmp_ans+=got[j];
                    continue;
                }
                while (mat.variable.find(make_parameter(x))!=mat.variable.end())
                {
                    x++;
                }
                string tmp=make_parameter(x),tmp2;
                x++;
                tmp2=to_string(mat.val[i][j]);
                got[j]=tmp;
                if(abs(mat.val[i][j])!=1)
                {
                    tmp_ans+=tmp2;
                }
                tmp_ans+=tmp;
            }
        }
        string tmp;
        if(chk[id]!=0)
        {
            tmp=to_string(chk[id]);
        }
        got[id]=tmp+tmp_ans;
    }
    x=0;
    for(auto var : mat.variable)
    {
        mat.ans[var]=got[x];
        x++;
    }
    return mat;
}

int main()
{
    string s;
    vector<string> arr;
    vector<int> c;
    set<string> variables;
    while (getline(cin, s))
    {
        s = make_normal(s);
        arr.push_back(s);
        stringstream exp(s);
        string part;
        while (exp >> part)
        {
            if (part[0] != '.' && !(part[0] >= '0' && part[0] <= '9') && part[0] != '+' && part[0] != '-' && part[0] != '=')
            {
                variables.insert(part);
            }
        }
    }
    MATRIX mat(arr.size(),variables.size() + 1,variables);

    for (int i = 0; i < (int)arr.size(); i++)
    {
        stringstream exp(arr[i]);
        map<string, double> got;
        string part;
        double d = 1, sign = 1, ok = 0;
        while (exp >> part)
        {
            if (part[0] == '+')
            {
                sign = 1;
            }
            else if (part[0] == '-')
            {
                sign = -1;
            }
            else if (part[0]=='.' || (part[0] >= '0' && part[0] <= '9'))
            {
                d = stod(part);
                if (ok)
                {
                    c.push_back(d*sign);
                }
            }
            else if (part[0] != '.' && !(part[0] >= '0' && part[0] <= '9') && part[0] != '+' && part[0] != '-' && part[0] != '=')
            {
                got[part] = d * sign;
                d = 1;
                sign = 1;
            }
            else if (part[0] == '=')
            {
                ok = 1;
            }
        }
        int j = 0;
        for (auto x : variables)
        {
            mat.val[i][j] = got[x];
            j++;
        }
    }
    cout<<"The Augmented Matrix is :"<<endl;

    for(int i=0;i<(int)c.size();i++)
    {
        mat.val[i][mat.m-1]=c[i];
    }
    mat.print();

    mat=make_row_echelon(mat);
    cout<<"Augmented Matrix in Row Echelon form is :"<<endl;
    mat.print();

    mat=make_reduced_row_echelon(mat);
    cout<<"Augmented Matrix in Reduced Row Echelon form is : "<<endl;
    mat.print();

    mat=make_ans(mat);
    cout<<"Solution according to the vector is : "<<endl;
    if(mat.exist)
    {
        for(auto var : mat.variable)
        {
            cout<<var<<" = "<<mat.ans[var]<<endl;
        }
    }
    else
    {
        cout<<"There is no solution exists"<<endl;
    }
    return 0;
}
