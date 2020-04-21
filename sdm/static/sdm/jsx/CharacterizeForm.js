class CharacterizeForm extends React.Component {
    constructor(props) {
        super(props)
        if (this.props.currentPage !== 'characterize') {
            return null;
        }
    }

    render() {
        return (
            <div className='form-group'>
                <label htmlFor='sequence'>Enter DNA sequence</label>
                <textarea
                    className='form-control md-textarea'
                    id='sequence'
                    name='sequence'
                    type='text'
                    value={this.props.sequence}
                    onChange={this.props.handleChange}
                    autofocus
                />
            </div>
        );
    }
}
