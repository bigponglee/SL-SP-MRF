function [RFpulses, TR] = generate_RF_TR(num_of_TR)
% load FA, TR
tmp = load('fatr.mat');
RFpulses = tmp.fa;
TR = tmp.tr;
RFpulses = RFpulses(1:num_of_TR) * pi / 180;
TR = TR(1:num_of_TR);
end