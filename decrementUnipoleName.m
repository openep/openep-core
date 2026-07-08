function uni2name = decrementUnipoleName(uni1name)
%increment the unipole name by 1 to get the second unipole (this assumes
%that the second unipole is always the 'first + 1'
pattern = '(?<name>\S*\D)(?<number>\d*)';
result = regexpi(uni1name,pattern,'names');
uni2name = [result.name, num2str(str2double(result.number)-1)];
end