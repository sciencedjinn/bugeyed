function R=cut(R)

azimin=-180;
azimax=180;
elemin=-58.75;
elemax=90;
R1=R(:, 1);
R2=R(:, 2);
R1(R1<azimin | R1>azimax | R2<elemin | R2>elemax)=66666;
R2=R2(R1~=66666);
R1=R1(R1~=66666);
R=cat(2, R1, R2);