%posted here http://www.mathworks.com/matlabcentral/fileexchange/13916--simple--tool-for-estimating-the-number-of-clusters/content/valid_RandIndex.m
%by by Kaijun Wang
%10 Feb 2007 (Updated 08 Jul 2009)
%patched to make fast enough for huge cluster results by C.Gaiteri 2014

%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the 
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of 
% Classification 2:193-218
%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk
%Copyright (c) 2006-2007, Kaijun WANG
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:%
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution

%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

function ARI=adjustedrand(partitionA,partitionB)

ng1=max(partitionA);
ng2=max(partitionB);
%ctabmat=full(sparse(partitionA,partitionB,1,ng1,ng2)); slightly faster, but memory intense for large matrices
ctabmat=(sparse(partitionA,partitionB,1,ng1,ng2));  %SpeakEasy passes in sequentially labeled clusters, so this shouldn't be an issue 

n=sum(sum(ctabmat));
nis=sum(sum(ctabmat,2).^2);		%sum of squares of sums of rows
njs=sum(sum(ctabmat,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(ctabmat.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   ARI=0;			%avoid division by zero; if k=1, define Rand = 0
else
   ARI=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end
