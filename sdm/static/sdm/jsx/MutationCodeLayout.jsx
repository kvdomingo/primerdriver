class MutationCodeLayout extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        if (this.props.mutation_type === 'sub') {
            return (
                <div className='row row-cols-1 row-cols-md-3'>
                    <div className='col form-group'>
                        <label htmlFor='target'>Target</label>
                        <input
                            type='text'
                            className='form-control'
                            name='target'
                            id='target'
                            value={this.props.target}
                            onChange={this.props.genericNucleobaseHandler}
                            required={true}
                        />
                    </div>
                    <div className='col form-group'>
                        <label htmlFor='target'>Position</label>
                        <input
                            type='number'
                            min='0'
                            className='form-control'
                            name='position'
                            id='position'
                            value={this.props.position}
                            onChange={this.props.genericIntHandler}
                            required={true}
                        />
                    </div>
                    <div className='col form-group'>
                        <label htmlFor='target'>Replacement</label>
                        <input
                            type='text'
                            className='form-control'
                            name='replacement'
                            id='replacement'
                            min='0'
                            value={this.props.replacement}
                            onChange={this.props.genericNucleobaseHandler}
                            required={true}
                        />
                    </div>
                </div>
            );
        } else if (this.props.mutation_type === 'ins') {
            return (
                <div className='row row-cols-1 row-cols-md-2'>
                    <div className='col form-group'>
                        <label htmlFor='target'>Position</label>
                        <input
                            type='number'
                            min='0'
                            className='form-control'
                            name='position'
                            id='position'
                            value={this.props.position}
                            onChange={this.props.genericIntHandler}
                            required={true}
                        />
                    </div>
                    <div className='col form-group'>
                        <label htmlFor='target'>Replacement</label>
                        <input
                            type='text'
                            className='form-control'
                            name='replacement'
                            id='replacement'
                            min='0'
                            value={this.props.replacement}
                            onChange={this.props.genericNucleobaseHandler}
                            required={true}
                        />
                    </div>
                </div>
            );
        } else if (this.props.mutation_type === 'del') {
            return (
                <div className='row row-cols-1 row-cols-md-2'>
                    <div className='col form-group'>
                        <label htmlFor='target'>Target</label>
                        <input
                            type='text'
                            className='form-control'
                            name='target'
                            id='target'
                            value={this.props.target}
                            onChange={this.props.genericNucleobaseHandler}
                            required={true}
                        />
                    </div>
                    <div className='col form-group'>
                        <label htmlFor='target'>Position</label>
                        <input
                            type='number'
                            min='0'
                            className='form-control'
                            name='position'
                            id='position'
                            value={this.props.position}
                            onChange={this.props.genericIntHandler}
                            required={true}
                        />
                    </div>
                </div>
            );
        } else {
            return(null)
        }
    }
}
