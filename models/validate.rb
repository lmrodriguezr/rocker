#!/usr/bin/env ruby

require 'json'

if ARGV.empty?
  $stderr.puts "Usage:"
  $stderr.puts "    #{__FILE__} list.json"
  exit(1)
end

list = JSON.parse(File.read(ARGV.first), symbolize_names: true)
raise "No families in the list" unless list[:families]

list[:families].each do |gene, family|
  # Evaluate each family
  [:product, :collections]. each do |k|
    warn "Family #{gene}: undefined #{k}" unless family[k]
  end
  warn "Family #{gene}: folder doesn't exist" unless
        Dir.exist? File.expand_path("../#{gene}", __FILE__)
  unk = family.keys - [:product, :collections]
  warn "Family #{gene}: Unknown keys: #{unk}" unless unk.empty?
  next unless family[:collections]

  family[:collections].each do |name, coll|
    # Evaluate each collection
    coll_name = "Collection #{gene}/#{name}"
    [:author, :update, :positive, :models].each do |k|
      warn "#{coll_name}: undefined #{k}" unless coll[k]
    end
    unk = coll.keys - [:description, :author, :update, :comments,
          :search, :positive, :negative, :models]
    warn "#{coll_name}: Unknown keys: #{unk}" unless unk.empty?
    %w[positive.txt negative.txt ref.fasta].each do |ext|
      next if ext=="negative.txt" and not coll[:negative]
      f = File.expand_path("../#{gene}/#{name}.#{ext}", __FILE__)
      warn "#{coll_name}: #{f} doesn't exist" unless File.exist? f
    end
    next unless coll[:models]

    coll[:models].each do |model|
      # Evaluate each model
      model_name = "Model #{gene}/#{name}/#{model[:length]}"
      warn "#{model_name}: undefined length" unless model[:length]
      [:stats, :'50bit_stats'].each do |k|
        if model[k] and model[k].size!=3
          warn "#{model_name}: #{k} is not a triple"
        end
      end
      unk = model.keys - [:length, :stats, :'50bit_stats']
      warn "#{model_name}: Unknown keys: #{unk}" unless unk.empty?
      %w[rocker jpg].each do |ext|
        f = File.expand_path("../#{gene}/#{name}.#{model[:length]}.#{ext}",
              __FILE__)
        warn "#{model_name}: #{f} doesn't exist" unless File.exist? f
      end
    end
  end
end

