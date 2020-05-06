import React from 'react';
import 'mdbreact';


export default function LoadingScreen() {
	return (
		<div className='container mb-5 text-center'>
			<div className='spinner-grow mb-2' role='status'></div>
			<p>Processing...</p>
		</div>
	);
}
