#!/usr/bin/perl -w
use strict;
use Module::Build;

my $build = Module::Build->new(
  module_name     => 'Math::MPath',
  license         => 'perl',
  requires => {
    'Math::BigFloat' => 0,
  },
  'build_requires' => {
    'Module::Build' => '0.3601',
    'Test::More' => '0',
  },
  'configure_requires' => {
    'Module::Build' => '0.38',
  },
);

$build->create_build_script;