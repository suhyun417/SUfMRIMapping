% remove scans in the old way:

rem=[1 1;16 17;19 21;26 27];

% fill up the removed scans vector
remove=zeros(1,numscans); k=0;
for i = 1:size(rem,1)
    for j = rem(i,1):rem(i,2)
        k=k+1;
        remove(k) = j;
    end
end
remove = remove(1:k);