function polyout = polybuffer_wrapper(list_of_inputs)
%trying out how to use polybuffer to create a buffer zone arund the lines
%list_of_inputs = {'lines_test';'lines_test_2'};
%what we have is a list of matrices, A, B, C, D.  The list is called
%Alphabet
%for each loop in the for loop, i need one of them
%function output = name of matlab script(Alphabet)

if ~exist('list_of_inputs','var')
    list_of_inputs = {'lines_test'};
end
polyout = cell(1,length(list_of_inputs));
for i = 1:length(list_of_inputs)
    load([list_of_inputs{i} '.mat'],list_of_inputs{i});
    polyout{i}= polybuffer(eval(list_of_inputs{i}),'points',0.03);
end

moo

%
    