package com.sanger.additionalBamStatistics;

import com.beust.jcommander.*;

public class PositiveInteger implements IParameterValidator {
	public void validate(String name, String value) throws ParameterException {
		int n = Integer.parseInt(value);
		if(n < 0)
			throw new ParameterException("Parameter " + name + " should be positive. Found " + value);
	}
}
