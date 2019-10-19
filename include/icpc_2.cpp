#include "bits/stdc++.h"
using namespace std;
int main()
{
	int t;
	cin>>t;
	while(t--)
	{
		int n,k;
		cin>>n>>k;
		int a[n];
		for(int i=0; i<n; i++)
			cin>>a[i];
		int curr_sum=a[0];
		int count = 0;
		for(int i=0; i<n-1; i++)
		{

			if(curr_sum > 0)
			{
				if(a[i+1] >= 0)
				{
					curr_sum = 0;
					count++;
				}
				else 
				{
					if(curr_sum+a[i+1]>0)
					{
						curr_sum+=a[i+1];
					}
					else
					{
						curr_sum = 0;
						count++;
					}
				}
			}
			else
			{

			}
		}
	}
}