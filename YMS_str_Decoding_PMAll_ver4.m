%% FASTQ 읽어오기
clear,clc
% 1. DNAdata Stores Apple
% 2. J_W Likes Orange
% 3. F_C Likes Grape
% 4. J_V_N Has Apple
% 5. DNAdata Stores Apple Orange
% 6. DNAdata Stores Parallelly [Apple, Orange, Grape]

filenumber = input("file 번호?: ",'s');
DNA1 = input("첫번째 단어: ",'s');
DNA2 = input("두번째 단어: ",'s');
DNA3 = input("세번째 단어: ",'s');
DNA4 = input("네번째 단어(없으면 그냥 엔터): ",'s');
PrimerEquality = input("프라이머 일치하는 개수?: ");
DNAinf = [{DNA1} DNA2 DNA3 DNA4];
[~,S1,Q1] = fastqread(filenumber+"_1.fastq");
[~,S2,Q2] = fastqread(filenumber+"_2.fastq");
%% file sequence
Primers1 = cell(1,2);
Primers2 = cell(1,2);
Primers3 = cell(1,2);
Primers4 = cell(1,2);

if sum(ismember(DNAinf,'DNAdata')) > 0
    eval("Primers"+find(ismember(DNAinf,'DNAdata'))+"{1} = 'TTCGTTCGTCGTTGATTGGT';")
    eval("Primers"+find(ismember(DNAinf,'DNAdata'))+"{2} = 'AAACGGAGCCATGAGTTTGT';")
end
if sum(ismember(DNAinf,'Stores')) > 0
    eval("Primers"+find(ismember(DNAinf,'Stores'))+"{1} = 'TCCTCAGCCGATGAAATTCC';")
    eval("Primers"+find(ismember(DNAinf,'Stores'))+"{2} = 'TGTACCATCCGTTTGACTGG';")
end
if sum(ismember(DNAinf,'Apple')) > 0
    eval("Primers"+find(ismember(DNAinf,'Apple'))+"{1} = 'CTGTCCATAGCCTTGTTCGT';")
    eval("Primers"+find(ismember(DNAinf,'Apple'))+"{2} = 'GCGGAAACGTAGTGAAGGTA';")
end
if sum(ismember(DNAinf,'J_W')) > 0
    eval("Primers"+find(ismember(DNAinf,'J_W'))+"{1} = 'GTCCAGGCAAAGATCCAGTT';")
    eval("Primers"+find(ismember(DNAinf,'J_W'))+"{2} = 'ACCACCGTTAGGCTAAAGTG';")
end
if sum(ismember(DNAinf,'Likes')) > 0
    eval("Primers"+find(ismember(DNAinf,'Likes'))+"{1} = 'TACCGCATCCTTATTCGAGC';")
    eval("Primers"+find(ismember(DNAinf,'Likes'))+"{2} = 'TCTGGTGCAAGCCAATGAAA';")
end
if sum(ismember(DNAinf,'Orange')) > 0
    eval("Primers"+find(ismember(DNAinf,'Orange'))+"{1} = 'TGTATTTCCTTCGGTGCTCC';")
    eval("Primers"+find(ismember(DNAinf,'Orange'))+"{2} = 'TTTCGACAACGGTCTGGTTT';")
end
if sum(ismember(DNAinf,'F_C')) > 0
    eval("Primers"+find(ismember(DNAinf,'F_C'))+"{1} = 'ATCCTGCAAACGCATTTCCT';")
    eval("Primers"+find(ismember(DNAinf,'F_C'))+"{2} = 'ATGCCTTTCCGAAGTTTCCA';")
end
if sum(ismember(DNAinf,'Grape')) > 0
    eval("Primers"+find(ismember(DNAinf,'Grape'))+"{1} = 'AGCCTTGTGTCCATCAATCC';")
    eval("Primers"+find(ismember(DNAinf,'Grape'))+"{2} = 'TGCGCTATGGTTTGGCTAAT';")
end
if sum(ismember(DNAinf,'J_V_N')) > 0
    eval("Primers"+find(ismember(DNAinf,'J_V_N'))+"{1} = 'TAGCCTCCAGAATGAAACGG';")
    eval("Primers"+find(ismember(DNAinf,'J_V_N'))+"{2} = 'TTCAAGCCAAACCGTGTGTA';")
end
if sum(ismember(DNAinf,'Has')) > 0
    eval("Primers"+find(ismember(DNAinf,'Has'))+"{1} = 'GAAGAGTTTAGCCACCTGGT';")
    eval("Primers"+find(ismember(DNAinf,'Has'))+"{2} = 'AAGGCCAATTCGCGGTTATT';")
end
if sum(ismember(DNAinf,'Parallelly')) > 0
    eval("Primers"+find(ismember(DNAinf,'Parallelly'))+"{1} = 'CAAGATTGTGGACGATTGGC';")
    eval("Primers"+find(ismember(DNAinf,'Parallelly'))+"{2} = 'TGCAATGTTTCCGTCGGTTT';")
end

%% file 1
TS1 = [S1 S2]; % 정방향 역방향 합치는것
if isempty(DNA4)
    Reverse = Primers3{2};
else
    Reverse = Primers4{2};
end
A_F = strfind(Reverse,'A');
T_F = strfind(Reverse,'T');
G_F = strfind(Reverse,'G');
C_F = strfind(Reverse,'C');
S_comp = Reverse;
S_comp(A_F) = 'T';
S_comp(T_F) = 'A';
S_comp(G_F) = 'C';
S_comp(C_F) = 'G';
rev = flip(S_comp);

S = cell(1,length(TS1));
for a = 1 : length(TS1)
    if sum(TS1{a}(1:length(rev)) == rev) >= length(rev)-2
        Reverse = TS1{a};
        A_F = strfind(Reverse,'A');
        T_F = strfind(Reverse,'T');
        G_F = strfind(Reverse,'G');
        C_F = strfind(Reverse,'C');
        S_comp = Reverse;
        S_comp(A_F) = 'T';
        S_comp(T_F) = 'A';
        S_comp(G_F) = 'C';
        S_comp(C_F) = 'G';
        S{a} = flip(S_comp);
    else
        S{a} = TS1{a};
    end
end

T = cell(length(S),1); % 빈 셀 생성
for a = 1 : length(S)
    if sum(S{a}(1:20) == Primers1{1}) >= PrimerEquality % 단어 F primer 비교, 18개 이상 일치하면 저장
        T{a} = S{a};
    end
end
Tempty = ~cellfun('isempty',T);
TT = T(Tempty); % 프라이머 분류 시퀀스
%%
R1 = cell(length(TT),4); % 빈 셀 생성 (세 줄)
PMidx = zeros(length(TT),4);
for b = 1 : length(TT)
    for c = 40 : length(TT{1})
        if sum(TT{b}(c-19:c) == Primers1{2}) >= PrimerEquality % 단어 R primer 비교
            R1{b,1} = TT{b}(1:c); % 첫번째 셀에 저장
            PMidx(b,1) = 1;
            break
        end
    end
    for d = 60 : length(TT{1})
        if sum(TT{b}(d-19:d) == Primers2{1}) >= PrimerEquality % 단어 F primer 비교
            for e = 80 : length(TT{1})
                if sum(TT{b}(e-19:e) == Primers2{2}) >= PrimerEquality % 단어 R primer 비교
                    R1{b,2} = TT{b}(d-19:e); % 두번째 셀에 저장
                    PMidx(b,2) = 1;
                break
                end
            end
            break
        end
    end
    for f = 100 : length(TT{1})
        if sum(TT{b}(f-19:f) == Primers3{1}) >= PrimerEquality % 단어 비교
            for g = 120 : length(TT{1})
                if sum(TT{b}(g-19:g) == Primers3{2}) >= PrimerEquality % 단어 비교
                    R1{b,3} = TT{b}(f-19:g); % 세번째 셀에 저장
                    PMidx(b,3) = 1;
                break
                end
            end
            break
        end
    end
    if ~isempty(DNA4)
        for h = 120 : length(TT{1})
            if sum(TT{b}(h-19:h) == Primers4{1}) >= PrimerEquality % 단어 비교
                for i = 100 : length(TT{1})
                    if sum(TT{b}(i-19:i) == Primers4{2}) >= PrimerEquality % 단어 비교
                        R1{b,4} = TT{b}(h-19:i); % 세번째 셀에 저장
                        PMidx(b,4) = 1;
                        break
                    end
                end
            end
            break
        end
    end
end

str1 = R1(:,1);
str2 = R1(:,2);
str3 = R1(:,3);
str4 = R1(:,4);
%% Primer All Matching
DNAdata = 'TTCGTTCGTCGTTGATTGGTGCGCTCACGGCAGTCGACAATCGATTAGTAACAAACGGAGCCATGAGTTTGT';
Stores = 'TCCTCAGCCGATGAAATTCCAATTTGCACTAGTGGTCCCCGGACCTTGTACCATCCGTTTGACTGG';
Apple = 'CTGTCCATAGCCTTGTTCGTCGGCGATGGTCACTCCATTATTCGCGCGGAAACGTAGTGAAGGTA';
J_W = 'GTCCAGGCAAAGATCCAGTTAAGTGGTCTATAGTTGTTTAACGACAATTGTTTCTATACCCCTGACACCACCGTTAGGCTAAAGTG';
Likes = 'TACCGCATCCTTATTCGAGCAAACCAGCAGCGAATTCCCCGGATACTCTGGTGCAAGCCAATGAAA';
Orange = 'TGTATTTCCTTCGGTGCTCCTTGCCATCGACTCTTCATAGTTACTCACTTTCGACAACGGTCTGGTTT';
F_C = 'ATCCTGCAAACGCATTTCCTACGCCATCGACTCTTGATGCTTGTTCGAAGTTGCTTGATCCTCTCGGGAATGCCTTTCCGAAGTTTCCA';
Grape = 'AGCCTTGTGTCCATCAATCCTCGCCATCGACGATCCATTAGATTGTGCGCTATGGTTTGGCTAAT';
J_V_N = 'TAGCCTCCAGAATGAAACGGAAGTTCTGAACTCTACCTTCTATACTCACCATAATCTATGGTATACTCTAGTGTTGTTCAAGCCAAACCGTGTGTA';
Has = 'GAAGAGTTTAGCCACCTGGTAACAGAGATTCCTGGTAAGGCCAATTCGCGGTTATT';
Parallelly = 'CAAGATTGTGGACGATTGGCAAAAATTAGTAGAAAAACCGAACTCGGAATTCCCCCCCCCGTTAGTTTGCAATGTTTCCGTCGGTTT';

PMsort = zeros(length(R1),1);
for a = 1 : length(R1)
    if ~isempty(str1(a))&&length(str1{a})==length(eval(DNA1))
        if ~isempty(str2(a))&&length(str2{a})==length(eval(DNA2))
            if ~isempty(str3(a))&&length(str3{a})==length(eval(DNA3))
                PMsort(a) = 1;
            end
        end
    end
end
%%
PMAll = R1(find(PMsort),:);

PMstr1 = PMAll(:,1);
PMstr2 = PMAll(:,2);
PMstr3 = PMAll(:,3);
PMstr4 = PMAll(:,4);

[PCodondata1,Pdata1] = strmode(PMstr1);
[PCodondata2,Pdata2] = strmode(PMstr2);
[PCodondata3,Pdata3] = strmode(PMstr3);
Pdatacell = [{Pdata1} {Pdata2} {Pdata3}];
if ~isempty(DNA4)
    [PCodondata4,Pdata4] = strmode(PMstr4);
    Pdatacell = [{Pdata1} {Pdata2} {Pdata3} {Pdata4}];
end

writematrix(PCodondata1,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"A2:D"+length(eval(DNA1))+1)
writematrix(PCodondata2,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"F2:I"+length(eval(DNA2))+1)
writematrix(PCodondata3,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"K2:N"+length(eval(DNA3))+1)

writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"A1:D1")
writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"F1:I1")
writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"K1:N1")
if ~isempty(DNA4)
    writematrix(PCodondata4,[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"P2:S"+length(eval(DNA4))+1)
    writematrix(["A" "T" "G" "C"],[DNA1 DNA2 DNA3 DNA4]+" PM All.xls",'Range',"P1:S1")
end
%% file Mode & Perfect Matching
[Codondata1,data1] = strmode(str1);
[Codondata2,data2] = strmode(str2);
[Codondata3,data3] = strmode(str3);
datacell = [{data1} {data2} {data3}];
DNA1+" 최빈값: "+ data1
DNA2+" 최빈값: "+ data2
DNA3+" 최빈값: "+ data3
csvwrite(DNA1+"_Codondata.csv",Codondata1)
csvwrite(DNA2+"_Codondata.csv",Codondata2)
csvwrite(DNA3+"_Codondata.csv",Codondata3)
if ~isempty(DNA4)
    [Codondata4,data4] = strmode(str4);
    datacell = [{data1} {data2} {data3} {data4}];
    DNA4+" 최빈값: "+ data4
    csvwrite(DNA4+"_Codondata.csv",Codondata4)
end

strmatching(DNAinf,datacell,'DNAdata')
strmatching(DNAinf,datacell,'Stores')
strmatching(DNAinf,datacell,'Apple')
strmatching(DNAinf,datacell,'J_W')
strmatching(DNAinf,datacell,'Likes')
strmatching(DNAinf,datacell,'Orange')
strmatching(DNAinf,datacell,'F_C')
strmatching(DNAinf,datacell,'Grape')
strmatching(DNAinf,datacell,'J_V_N')
strmatching(DNAinf,datacell,'Has')
strmatching(DNAinf,datacell,'Parallelly')