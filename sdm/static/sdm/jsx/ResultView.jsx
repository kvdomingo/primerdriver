class ResultView extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        if (typeof(this.props.results) === 'object') {
            if (this.props.mode === 'CHAR') {
                table_content = []
                for (key in this.props.results) {
                    table_content.push(
                        <React.Fragment key={key}>
                            <tr>
                                <th scope='row'>{key}</th>
                                <td>{this.props.results[key]['1']}</td>
                            </tr>
                        </React.Fragment>
                    )
                }
                return (
                    <div className='table-responsive text-nowrap'>
                        <a
                            className='btn btn-blue-grey mb-4 mr-3'
                            id='back'
                            onClick={function(e) {this.props.responseCatcher(e, null); this.props.changeView(e, 0)}.bind(this)}
                        >
                            <i className='fas fa-arrow-left mr-3'></i>main menu
                        </a>
                        <table className='table table-sm table-bordered table-hover'>
                            <tbody>{table_content}</tbody>
                        </table>
                    </div>
                );
            } else {
                return ((() => {
                    table_content = [
                        <React.Fragment>
                            <a
                                className='btn btn-blue-grey mb-4 mr-3'
                                id='back'
                                onClick={function(e) {this.props.responseCatcher(e, null); this.props.changeView(e, 0)}.bind(this)}
                            >
                                <i className='fas fa-arrow-left mr-3'></i>main menu
                            </a>
                            <h2 className='text-md-center mx-md-2 h2-responsive d-md-inline'>{Object.keys(this.props.results).length} results</h2>
                        </React.Fragment>
                    ];
                    for (let i = 0, n = Object.keys(this.props.results).length; i < n; i++) {
                        table_content.push(
                            <React.Fragment key={i}>
                                <div className='table-responsive text-nowrap'>
                                    <table className='table table-sm table-bordered table-hover'>
                                        <thead>
                                            <tr>
                                                <th scope='col'></th>
                                                <th scope='col'>{`Primer ${i+1}`}</th>
                                            </tr>
                                        </thead>
                                        <tbody>{(() => {
                                            row_content = [];
                                            for (let key in this.props.results[(i+1).toString()]) {
                                                row_content.push(
                                                    <React.Fragment key={key}>
                                                        <tr>
                                                            <th scope='row'>{key}</th>
                                                            <td>{this.props.results[(i+1).toString()][key]}</td>
                                                        </tr>
                                                    </React.Fragment>
                                                );
                                            }
                                            return row_content;
                                        })()}
                                        </tbody>
                                    </table>
                                </div>
                            </React.Fragment>
                        );
                    }
                    return table_content;
                })());
            }
        } else {
            return (
                <div className='text-center'>
                    <a
                        className='btn btn-blue-grey btn-rounded mb-4 mr-3'
                        href='/'
                        id='back'
                        onClick={(e) => this.props.changeView(e, 0)}
                    >
                        <i className='fas fa-arrow-left mr-3'></i>main menu
                    </a>
                    <p>An error occurred. Please try again later, or <a href='https://github.com/kvdomingo/primerdriver/issues' target='_blank'>report an issue</a>.</p>
                </div>
            );
        }
    }
}
